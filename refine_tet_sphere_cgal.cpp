#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;


typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;
// // Criteria
// typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// Sizing field
struct Spherical_sizing_field
{
  typedef ::FT FT;
  typedef Point Point_3;
  typedef Mesh_domain::Index Index;
  
  FT operator()(const Point_3& p, const int, const Index&) const
  {
    FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
    return CGAL::abs( CGAL::sqrt(sq_d_to_origin)-0.5 ) / 100. + 0.025; 
  }
};

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input("data/elephant.off");
  input >> polyhedron;
   
  // Create domain
  Mesh_domain domain(polyhedron);
  
  Spherical_sizing_field size;
  // Mesh criteria (no cell_size set)
  Facet_criteria facet_criteria(30, 0.1, 0.025); // angle, size, approximation
  Cell_criteria cell_criteria(2, size); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);
  // Mesh_criteria criteria(facet_angle=25, facet_size=0.15, facet_distance=0.008,
  //                          cell_radius_edge=3);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  // Output
  std::ofstream medit_file("out_1.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();

  // // Set tetrahedron size (keep cell_radius_edge), ignore facets
  //  Mesh_criteria new_criteria(cell_radius_edge=3, cell_size=0.03);
  // 
  //  // Mesh refinement
  //  CGAL::refine_mesh_3(c3t3, domain, new_criteria);
  // 
  //  // Output
  //  medit_file.open("out_2.mesh");
  //  c3t3.output_to_medit(medit_file);

  return 0;
}
