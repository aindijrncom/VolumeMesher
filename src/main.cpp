#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "BSP.h"

/// <summary>
/// Read a triangle mesh from an OFF file.
/// Robust parser: skips comments (#), empty lines;
/// accepts "OFF", "OFF BINARY" (text mode only);
/// handles missing edge count, vertex colors, face colors.
/// </summary>
void read_OFF_file(const char* filename,
    double** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose) {

    FILE* file = fopen(filename, "r");
    if (file == NULL)
        ip_error("Cannot open input file.\n");

    char line[4096];
    int found_off = 0;

    // Phase 1: find "OFF" header, skipping comments and blank lines
    while (fgets(line, sizeof(line), file)) {
        // Trim leading whitespace
        char* p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '\n' || *p == '\r') continue; // blank line
        if (*p == '#') continue; // comment
        // Skip UTF-8 BOM if present (Windows/Meshlab exports)
        if ((unsigned char)p[0] == 0xEF &&
            (unsigned char)p[1] == 0xBB &&
            (unsigned char)p[2] == 0xBF) {
            p += 3;
            // If BOM is the only content on this line, skip it
            while (*p && isspace((unsigned char)*p)) p++;
            if (*p == '\0' || *p == '\n' || *p == '\r') continue;
        }
        // Check if header ends with OFF (handles OFF, NOFF, COFF, STOFF, etc.)
        size_t hlen = strcspn(p, " \t\r\n");
        if (hlen >= 3 && strncmp(p + hlen - 3, "OFF", 3) == 0) { found_off = 1; break; }
        ip_error("Invalid input file format: expected OFF header.\n");
    }
    if (!found_off) ip_error("Invalid input file format: no OFF header.\n");

    // Phase 2: read vertex/face count (skip comments/blank lines)
    while (fgets(line, sizeof(line), file)) {
        char* p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '\n' || *p == '\r') continue;
        if (*p == '#') continue;
        // Try 3 numbers (nv nf ne), fall back to 2 (nv nf)
        unsigned int _ne;
        int r = sscanf(p, "%u %u %u", npts, ntri, &_ne);
        if (r == 2 || r == 3) break;
        ip_error("Invalid input file data: expected vertex/face count.\n");
    }

    if (verbose)
        printf("file %s contains %u vertices and %u triangles\n",
               filename, *npts, *ntri);

    // Allocate storage
    *vertices_p = (double*)malloc(sizeof(double) * 3 * (*npts));
    *tri_vertices_p = (uint32_t*)malloc(sizeof(uint32_t) * 3 * (*ntri));

    // Phase 3: read vertex coordinates (skip comments)
    uint32_t v_read = 0;
    while (v_read < *npts && fgets(line, sizeof(line), file)) {
        char* p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '\n' || *p == '\r') continue;
        if (*p == '#') continue;
        double* v = (*vertices_p) + v_read * 3;
        if (sscanf(p, "%lf %lf %lf", v, v + 1, v + 2) < 3)
            ip_error("Invalid vertex data in input file.\n");
        v_read++;
    }
    if (v_read != *npts)
        ip_error("Unexpected end of file while reading vertices.\n");

    // Phase 4: read face indices (skip comments)
    uint32_t f_read = 0;
    while (f_read < *ntri && fgets(line, sizeof(line), file)) {
        char* p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '\n' || *p == '\r') continue;
        if (*p == '#') continue;
        unsigned int nv;
        uint32_t* f = (*tri_vertices_p) + f_read * 3;
        int r = sscanf(p, "%u %u %u %u", &nv, f, f + 1, f + 2);
        if (r < 4) ip_error("Invalid face data in input file.\n");
        if (nv != 3) ip_error("Non-triangular faces not supported.\n");
        f_read++;
    }
    if (f_read != *ntri)
        ip_error("Unexpected end of file while reading faces.\n");

    fclose(file);
}

namespace
{
uint32_t read_uint32_le(const unsigned char* data)
{
    return static_cast<uint32_t>(data[0]) |
           (static_cast<uint32_t>(data[1]) << 8) |
           (static_cast<uint32_t>(data[2]) << 16) |
           (static_cast<uint32_t>(data[3]) << 24);
}

float read_float32_le(const unsigned char* data)
{
    static_assert(sizeof(float) == sizeof(uint32_t), "STL reader requires 32-bit float values.");
    const uint32_t bits = read_uint32_le(data);
    float value = 0.0f;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

bool is_binary_STL_file(const char* filename, uint32_t* triangle_count)
{
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file)
        ip_error("Cannot open input file.\n");

    const std::streamoff file_size = file.tellg();
    if (file_size < 84)
        return false;

    std::array<unsigned char, 4> count_bytes{};
    file.seekg(80, std::ios::beg);
    file.read(reinterpret_cast<char*>(count_bytes.data()), count_bytes.size());
    if (!file)
        return false;

    const uint32_t count = read_uint32_le(count_bytes.data());
    const uint64_t expected_size = 84ULL + 50ULL * count;
    if (expected_size > static_cast<uint64_t>(file_size))
        return false;

    *triangle_count = count;
    return true;
}

void read_binary_STL_coordinates(const char* filename, uint32_t triangle_count, std::vector<double>* coordinates)
{
    if (triangle_count > std::numeric_limits<uint32_t>::max() / 3U)
        ip_error("STL file contains too many triangles.\n");

    std::ifstream file(filename, std::ios::binary);
    if (!file)
        ip_error("Cannot open input file.\n");

    file.seekg(84, std::ios::beg);
    coordinates->reserve(static_cast<size_t>(triangle_count) * 9U);

    std::array<unsigned char, 50> record{};
    for (uint32_t triangle = 0; triangle < triangle_count; ++triangle)
    {
        file.read(reinterpret_cast<char*>(record.data()), record.size());
        if (!file)
            ip_error("Unexpected end of binary STL file.\n");

        for (uint32_t vertex = 0; vertex < 3; ++vertex)
        {
            for (uint32_t component = 0; component < 3; ++component)
            {
                const size_t offset = 12U + vertex * 12U + component * 4U;
                const float value = read_float32_le(record.data() + offset);
                if (!std::isfinite(value))
                    ip_error("Binary STL contains a non-finite vertex coordinate.\n");
                coordinates->push_back(static_cast<double>(value));
            }
        }
    }
}

void read_ascii_STL_coordinates(const char* filename, std::vector<double>* coordinates)
{
    std::ifstream file(filename);
    if (!file)
        ip_error("Cannot open input file.\n");

    std::string token;
    while (file >> token)
    {
        std::transform(token.begin(), token.end(), token.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        if (token != "vertex")
            continue;

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(file >> x >> y >> z))
            ip_error("Invalid vertex data in ASCII STL file.\n");
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            ip_error("ASCII STL contains a non-finite vertex coordinate.\n");

        coordinates->push_back(x);
        coordinates->push_back(y);
        coordinates->push_back(z);
    }

    if (coordinates->empty() || coordinates->size() % 9U != 0)
        ip_error("Invalid ASCII STL file: expected three vertices per triangle.\n");
}

void copy_STL_triangle_soup(const std::vector<double>& coordinates,
    double** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri)
{
    const size_t vertex_count = coordinates.size() / 3U;
    const size_t triangle_count = vertex_count / 3U;
    if (vertex_count > std::numeric_limits<uint32_t>::max() ||
        triangle_count > std::numeric_limits<uint32_t>::max())
        ip_error("STL file is too large for 32-bit mesh indices.\n");

    *npts = static_cast<uint32_t>(vertex_count);
    *ntri = static_cast<uint32_t>(triangle_count);
    *vertices_p = static_cast<double*>(malloc(coordinates.size() * sizeof(double)));
    *tri_vertices_p = static_cast<uint32_t*>(malloc(vertex_count * sizeof(uint32_t)));
    if (*vertices_p == NULL || *tri_vertices_p == NULL)
    {
        free(*vertices_p);
        free(*tri_vertices_p);
        *vertices_p = NULL;
        *tri_vertices_p = NULL;
        ip_error("Not enough memory to load STL file.\n");
    }

    std::copy(coordinates.begin(), coordinates.end(), *vertices_p);
    for (uint32_t vertex = 0; vertex < *npts; ++vertex)
        (*tri_vertices_p)[vertex] = vertex;
}

void read_STL_file(const char* filename,
    double** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose)
{
    uint32_t binary_triangle_count = 0;
    const bool binary = is_binary_STL_file(filename, &binary_triangle_count);

    std::vector<double> coordinates;
    if (binary)
        read_binary_STL_coordinates(filename, binary_triangle_count, &coordinates);
    else
        read_ascii_STL_coordinates(filename, &coordinates);

    copy_STL_triangle_soup(coordinates, vertices_p, npts, tri_vertices_p, ntri);

    if (verbose)
        printf("file %s contains %u STL triangles (%s)\n",
            filename, *ntri, binary ? "binary" : "ASCII");
}

void read_mesh_file(const char* filename,
    double** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose)
{
    std::string extension(filename);
    const size_t dot = extension.find_last_of('.');
    extension = (dot == std::string::npos) ? std::string() : extension.substr(dot);
    std::transform(extension.begin(), extension.end(), extension.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (extension == ".off")
        read_OFF_file(filename, vertices_p, npts, tri_vertices_p, ntri, verbose);
    else if (extension == ".stl")
        read_STL_file(filename, vertices_p, npts, tri_vertices_p, ntri, verbose);
    else
        ip_error("Unsupported input format. Use an OFF or STL file.\n");
}
}

/// <summary>
/// Main function
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3) {
        printf("\nUsage:\n"
               "  mesh_generator input.{off|stl}\n"
               "  mesh_generator inputA.{off|stl} inputB.{off|stl}\n\n"
               "With one input, tetrahedralizes the enclosed volume.\n"
               "With two inputs, tetrahedralizes their union.\n\n"
               "OFF and STL inputs can be mixed; ASCII and binary STL are supported.\n\n"
               "Outputs:\n"
               "  skin.off    - triangulated boundary surface\n"
               "  skin.stl    - high-precision ASCII STL boundary surface\n"
               "  skin_outer.off - outer boundary only (enclosed cavities excluded)\n"
               "  skin_outer.stl - high-precision ASCII STL outer boundary only\n"
               "  volume.vtu  - tetrahedral volume mesh (VTU format)\n"
               "  volume.inp  - Abaqus/CalculiX C3D4 tetrahedral mesh\n\n"
               "Examples:\n"
               "  mesh_generator model.stl\n"
               "  mesh_generator ant.off pig.stl\n");
        return EXIT_FAILURE;
    }

    const bool two_input = (argc == 3);
    const char* fileA_name = argv[1];
    const char* fileB_name = two_input ? argv[2] : NULL;
    const char bool_opcode = two_input ? 'U' : '0';

    if (two_input)
        printf("Loading %s and %s (union)\n", fileA_name, fileB_name);
    else
        printf("Loading %s (single model)\n", fileA_name);

    double* coords_A = NULL;
    double* coords_B = NULL;
    uint32_t ncoords_A = 0;
    uint32_t ncoords_B = 0;
    uint32_t* tri_idx_A = NULL;
    uint32_t* tri_idx_B = NULL;
    uint32_t ntriidx_A = 0;
    uint32_t ntriidx_B = 0;

    read_mesh_file(fileA_name, &coords_A, &ncoords_A, &tri_idx_A, &ntriidx_A, false);
    if (two_input)
        read_mesh_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, false);

    BSPcomplex* complex = makePolyhedralMesh(
        fileA_name, coords_A, ncoords_A, tri_idx_A, ntriidx_A,
        fileB_name, coords_B, ncoords_B, tri_idx_B, ntriidx_B,
        bool_opcode, true, false, false);

    printf("Saving skin_outer.off and skin_outer.stl ...\n");
    complex->saveOuterSkin("skin_outer.off", bool_opcode, true, "skin_outer.stl");

    printf("Saving volume.vtu and volume.inp ...\n");
    complex->saveVolumeFiles("volume.vtu", "volume.inp", bool_opcode);

    printf("Saving skin.off and skin.stl ...\n");
    complex->saveSkin("skin.off", bool_opcode, true, "skin.stl");

    printf("Done.\n");
    return 0;
}
