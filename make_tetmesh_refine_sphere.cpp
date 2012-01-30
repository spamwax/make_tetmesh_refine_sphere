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
    double ref_ratio, def_size;
    double R1, R2, S1, S2;
    
    double A[3], B[3];
	double v2[3];
	double AB, m1, m2;
	
    FT operator()(const Point_3& p, const int, const Index&) const
    {
        // FT sq_d = CGAL::squared_distance(p, Point(x,y,z));
        // sq_d = CGAL::sqrt(sq_d);
        double P[3];
		double v1[3];
        P[0] = p.x(); P[1] = p.y(); P[2] = p.z();
        v_make(A, P, 3, v1); // v1 = P-A
        // v_make(A, B, 3, v2); // v2 = B-A
        
        double k = v_dot(v1, v2, 3) / AB / AB;
        double pp[3] = {k*v2[0]+A[0], k*v2[1]+A[1], k*v2[2]+A[2]};

        double pp_dist = v_dist(pp, P, 3);
        double localr = m1 * k * AB + R2;

        double localsize = m2 * k * AB + S2;
        // Check if point P lies within the 'cone'
        if (pp_dist <= localr) {
            return localsize;
        }
        else {
            return S2; //fabs(pp_dist - localr)/localr 
        }
	}
};



using namespace CGAL::parameters;

int main(int argc, char *argv[]) {
    Polyhedron cutbrain;
	std::string ifn, cfn, outfn;
    double facet_angle=25, facet_size=3.0, facet_distance=0.5,
           cell_radius_edge=3,cell_size=3.0;
    double xa = 0, ya = 0, za = 0, r1 = 1., ref_ratio = 0.3;
    double xb = 0, yb = 0, zb = 0, r2 = 1.;
    double s1 = 1., s2 = 1.;
    
    bool defulatcriteria = false;

    if (argc == 1) {
        std::cout << " Enter the surface mesh filename: ";
        std::cin >> ifn;
		std::cout << " Enter the output mesh file name: ";
		std::cin >> outfn;
        defulatcriteria = true;
        std::cout << " (Using default settings for meshing parameters!)\n";
    }
    else if (argc==3) {
        ifn = argv[1];
		outfn = argv[2];
        defulatcriteria = true;
    }
    else if (argc == 4) {
        ifn = argv[1];
		outfn = argv[2];
        cfn = argv[3];
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
        cfs >> xa; // Point A (tumor side)
        cfs >> ya;
        cfs >> za;
        cfs >> xb; // Point B (brain side)
        cfs >> yb;
        cfs >> zb;
        cfs >> r1; // radius on tumor side
        cfs >> r2; // radius on brain side
        cfs >> s1; // tet size on tumor side 
        cfs >> s2; // tet size on brain side
        cfs >> ref_ratio;
    }
    
    std::ifstream ifs(ifn.c_str());
    if (!ifs) {
        std::cerr << " Could not open input surface mesh file!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "    cgal_mesh: reading the input polyhedron.\n";
    ifs >> cutbrain;
//    std::cout << "    cgal_mesh: done reading the input polyhedron.\n";
    
    if (!cutbrain.is_closed()) {
        std::cerr << " Input polyhedron is not closed!\n";
        exit(-1);
    }
    
    Mesh_domain domain(cutbrain);

//    std::cout << "    cgal_mesh: setting domain criteria.\n";
    
    Spherical_sizing_field2 size;
    
    size.ref_ratio = ref_ratio;
    size.B[0] = xb; size.B[1] = yb; size.B[2] = zb;
    size.A[0] = xa; size.A[1] = ya; size.A[2] = za;
    size.def_size = cell_size;
    size.R1 = r1; size.R2 = r2;
    size.S1 = s1; size.S2 = s2;
    v_make(size.A, size.B, 3, size.v2);
	size.AB = v_magn(size.v2, 3);
	size.m1 = (size.R1 - size.R2) / size.AB;
	size.m2 = (size.S1 - size.S2) / size.AB;
	
    std::cout << "facet_angle " << facet_angle << ' ' << std::endl;
    std::cout << "facet_size " << facet_size << ' ' << std::endl;
    std::cout << "facet_distance " << facet_distance << ' ' << std::endl;
    std::cout << "cell_radius_edge " << cell_radius_edge << ' ' << std::endl;
    std::cout << "cell_size " << cell_size << ' ' << std::endl;
    std::cout << "A: " << xa << ' ' << ya << ' ' << za << std::endl <<
                 "B: " << xb << ' ' << yb << ' ' << zb << std::endl <<
                 "R1: " << r1 << std::endl <<
                 "R2: " << r2 << std::endl <<
                 "S1: " << s1 << std::endl <<
                 "S2: " << s2 << std::endl << 
                 "ref_ratio: " << ref_ratio << std::endl;
                 
    // Mesh criteria (no cell_size set)
    Facet_criteria facet_criteria(facet_angle, facet_size, facet_distance); // angle, size, approximation
    Cell_criteria cell_criteria(cell_radius_edge, size); // radius-edge ratio, size
    Mesh_criteria criteria(facet_criteria, cell_criteria);
    // Mesh_criteria criteria(facet_angle, facet_size, facet_distance, cell_radius_edge, cell_size);

    std::cout << "    cgal_mesh: running make_mesh_3.\n";
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

    std::ofstream medit_file(outfn.c_str());
    c3t3.output_to_medit(medit_file);
    medit_file.close();
  
    return EXIT_SUCCESS;
    
}

