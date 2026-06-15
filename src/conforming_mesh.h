#ifndef _CONFORMING_MESH_
#define _CONFORMING_MESH_

#include "delaunay.h"

#define CONSTR_GROUP_T uint32_t
#define CONSTR_A 0
#define CONSTR_B 1

// half-edges struct
struct half_edge_t {
  uint32_t endpts[2]; // endpoints of a constraint(triangle) edge,
                      // be sure that endpts[0] < endpts[1].
  uint32_t tri_ind;   // id_number of the constraint: tri_vertices[i]/3, where
                      // tri_vertices[i]=endpts[0] or tri_vertices[i]=endpts[1].
};

// constraints are triangles
class constraints_t {
public:
  uint32_t* tri_vertices;
  uint32_t num_triangles;
  uint32_t* constr_group;
  uint32_t num_virtual_triangles; // virtual constraints are indexed
                                  // from num_triangles - num_virtual_triangles
                                  // to num_triangles-1.

  constraints_t() : tri_vertices(NULL), num_triangles(0), constr_group(NULL), num_virtual_triangles(0) {};
  ~constraints_t()
  {
      if (tri_vertices) free(tri_vertices);
      if (constr_group) free(constr_group);
  }
};

// Edge constraints: independent of triangle constraints_t.
// edge_vertices stores [v0,v1, v0,v1, ...] for each edge,
// pointing into mesh->vertices[] after both remaps are applied.
struct edge_constraints_t {
    uint32_t* edge_vertices;  // length = 2 * num_edges
    uint32_t num_edges;

    edge_constraints_t() : edge_vertices(NULL), num_edges(0) {}
    ~edge_constraints_t() { if (edge_vertices) free(edge_vertices); }
};

void fill_half_edges(const constraints_t* constraints, half_edge_t* half_edges);
void sort_half_edges(half_edge_t* half_edges, uint32_t num_half_edges);
uint32_t place_virtual_constraints(TetMesh* mesh, constraints_t* constraints,
                               half_edge_t* half_edges);
void insert_constraints( TetMesh* , constraints_t*, uint32_t*, uint32_t**,
                                                   uint32_t*, uint32_t**,
                                                   uint32_t*, uint32_t**,
                                                   uint32_t*, uint32_t**,
                                                   uint32_t*, uint32_t**  );

// --- MVP-0: edge segment tracing ---

// Trace a single segment <v_start, v_stop> through the tetrahedral mesh.
// Returns touched tets via touched_tets / num_touched_tets (caller frees).
// mark_TetIntersection must be pre-allocated (calloc) by the caller,
// length = mesh->tet_num, and will be used as 0/INTERSECTION visited marker.
void trace_segment_through_tets(
    TetMesh* mesh,
    uint32_t v_start,
    uint32_t v_stop,
    uint32_t* mark_TetIntersection,
    uint64_t** touched_tets,
    uint64_t* num_touched_tets
);

// Trace all edges in edge_constraints and print debug output.
// Does NOT modify mesh or constraint maps.
void debug_trace_edge_constraints(
    TetMesh* mesh,
    edge_constraints_t* edge_constraints
);

#endif
