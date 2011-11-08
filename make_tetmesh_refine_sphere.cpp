#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>

#include "vector.h"

// Computational Kernel Definitions necessary for robust meshing
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
//typedef CGAL::Polyhedron_3<K>      Polyhedron;
//typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron> Mesh_domain;

// Domain
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Polyhedron_3<K> Polyhedron;
// typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
typedef CGAL::Polyhedron_3<Geom_traits> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Geom_traits> Mesh_domain;


typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);

// Meshing typedefs
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

// Sizing field
struct Spherical_sizing_field
{
    typedef ::FT FT;
    typedef Point Point_3;
    typedef Mesh_domain::Index Index;
    double x, y, z;
    double r, ref_ratio, def_size;

    FT operator()(const Point_3& p, const int, const Index&) const
    {
        FT sq_d = CGAL::squared_distance(p, Point(x,y,z));
        sq_d = CGAL::sqrt(sq_d);
        if (sq_d <= r + 0.01)
			return ref_ratio * def_size;
        else
            return def_size;
    // FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
    // return CGAL::abs( CGAL::sqrt(sq_d_to_origin)-0.5 ) / 10. + 0.025; 
    }
};

struct Spherical_sizing_field2
{
    typedef ::FT FT;
    typedef Point Point_3;
    typedef Mesh_domain::Index Index;
//  double x, y, z;
    double r, ref_ratio, def_size;

    double A[3], B[3];

    FT operator()(const Point_3& p, const int, const Index&) const
    {
      // FT sq_d = CGAL::squared_distance(p, Point(x,y,z));
      // sq_d = CGAL::sqrt(sq_d);
      double P[3];
      double v1[3], v2[3];
      P[0] = p.x(); P[1] = p.y(); P[2] = p.z();
      v_make(A, P, 3, v1);
      v_make(A, B, 3, v2);
      double k = v_dot(v1, v2, 3) / v_magn(v2, 3) / v_magn(v2, 3);

      bool within = false;
      if (k>=0 && k<=1)
            within = true;
        
      if (within)
          return k * def_size;
      else
          return def_size;
    // FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
    // return CGAL::abs( CGAL::sqrt(sq_d_to_origin)-0.5 ) / 10. + 0.025; 
	}
};



using namespace CGAL::parameters;

int main(int argc, char *argv[]) {
    Polyhedron cutbrain;
    std::string ifn, cfn;
    double facet_angle=25, facet_size=3.0, facet_distance=0.5,
           cell_radius_edge=3,cell_size=3.0;
    double x = 0, y = 0, z = 0, r = 1., ref_ratio = 0.3;
    
    bool defulatcriteria = false;

    if (argc == 1) {
        std::cout << " Enter the surface mesh filename: ";
        std::cin >> ifn;
        defulatcriteria = true;
        std::cout << " (Using default settings for meshing parameters!)\n";
    }
    else if (argc==2) {
        ifn = argv[1];
        defulatcriteria = true;
    }
    else if (argc == 3) {
        ifn = argv[1];
        cfn = argv[2];
    }

    if (!defulatcriteria) {
        std::cout << " Reading criteria from file." << std::endl;
        std::ifstream cfs(cfn.c_str());
        if (!cfs) {
            std::cerr << " Can not read mesh criteria file!\n";
            exit(-1);
        }
        cfs >> facet_angle;
        cfs >> facet_size;
        cfs >> facet_distance;
        cfs >> cell_radius_edge;
        cfs >> cell_size;
        cfs >> x;
        cfs >> y;
        cfs >> z;
        cfs >> r;
        cfs >> ref_ratio;
    }
    
    std::ifstream ifs(ifn.c_str());
    if (!ifs) {
        std::cerr << " Could not open input surface mesh file!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "    cgal_mesh: reading the input polyhedron.\n";
    ifs >> cutbrain;
    std::cout << "    cgal_mesh: done reading the input polyhedron.\n";
    
    if (!cutbrain.is_closed()) {
        std::cerr << " Input polyhedron is not closed!\n";
        exit(-1);
    }
    
    Mesh_domain domain(cutbrain);

    std::cout << "    cgal_mesh: setting domain criteria.\n";
    
    Spherical_sizing_field2 size;
    
    size.ref_ratio = ref_ratio;
    size.r = r;
    size.B[0] = x; size.B[1] = y; size.B[2] = z;
    size.A[0] = x+50; size.A[1] = y+50; size.A[2] = z+50;
    size.def_size = cell_size;
    
    // Mesh criteria (no cell_size set)
    Facet_criteria facet_criteria(facet_angle, facet_size, facet_distance); // angle, size, approximation
    Cell_criteria cell_criteria(cell_radius_edge, size); // radius-edge ratio, size
    Mesh_criteria criteria(facet_criteria, cell_criteria);
    // Mesh_criteria criteria(facet_angle, facet_size, facet_distance, cell_radius_edge, cell_size);

    std::cout << "    cgal_mesh: running make_mesh_3.\n";
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

    std::ofstream medit_file("out.mesh");
    c3t3.output_to_medit(medit_file);
    medit_file.close();
    std::cout << facet_angle << ' ' << std::endl;
    std::cout << facet_size << ' ' << std::endl;
    std::cout << facet_distance << ' ' << std::endl;
    std::cout << cell_radius_edge << ' ' << std::endl;
    std::cout << cell_size << ' ' << std::endl;
    std::cout << x << ' ' << y << ' ' << z << ' ' << r << std::endl;
    std::cout << ref_ratio << std::endl;
    return EXIT_SUCCESS;
    
}

