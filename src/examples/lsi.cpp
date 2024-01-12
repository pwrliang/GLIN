#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/LinearRing.h>
#include "omp.h"
#include "../stopwatch.h"
#include "alex.h"
#include <map>
#include "../../glin/glin.h"
#include "../../glin/piecewise.h"
#include "read_wkt.h"
//#include "../glin_benchmark/utils.h"
#include <geos/index/strtree/SimpleSTRtree.h>
#include <geos/index/strtree/GeometryItemDistance.h>
#include <geos/index/ItemVisitor.h>
#include <geos/geom/Envelope.h>
#include <geos/index/quadtree/Quadtree.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> point_type;
typedef boost::geometry::model::linestring<point_type> linestring_type;
typedef boost::geometry::model::polygon<point_type> polygon_type;
typedef boost::geometry::model::segment<point_type> segment_type;
typedef boost::geometry::model::multi_linestring<linestring_type> multi_lingstring_type;
typedef boost::geometry::model::multi_polygon<polygon_type> mpolygon_type;
typedef boost::variant<point_type, linestring_type, polygon_type, segment_type, multi_lingstring_type, mpolygon_type> geometryVariant;

geos::geom::Geometry *CreateQuery(double minX,
                                  double minY,
                                  double maxX,
                                  double maxY,
                                  geos::geom::GeometryFactory::Ptr &global_factory) {
    geos::geom::CoordinateArraySequence *cl6 = new geos::geom::CoordinateArraySequence();
    cl6->add(geos::geom::Coordinate(minX, minY));
    cl6->add(geos::geom::Coordinate(minX, maxY));
    cl6->add(geos::geom::Coordinate(maxX, maxY));
    cl6->add(geos::geom::Coordinate(maxX, minY));
    cl6->add(geos::geom::Coordinate(minX, minY));
    geos::geom::LinearRing *lr6 = global_factory->createLinearRing(cl6);
    geos::geom::Polygon *poly6 = global_factory->createPolygon(lr6, NULL);
    return poly6;
}


int main(int arc, char **argv) {
    alex::Glin<double, geos::geom::Geometry *> index;

    std::string map1(argv[1]);
    std::string map2(argv[2]);

    auto line_strings = ReadWKT(map1);
    auto queries = ReadWKT(map2);

    std::vector<std::unique_ptr<geos::geom::Geometry>> base_vec;
    std::vector<geos::geom::Geometry *> base_ptrs;
    std::vector<std::tuple<double, double, double, double>> pieces;

    double error_bound = 10000;
    std::unique_ptr<geos::geom::PrecisionModel> pm(new geos::geom::PrecisionModel());
    geos::geom::GeometryFactory::Ptr global_factory = geos::geom::GeometryFactory::create(pm.get(), -1);
    size_t total_segs = 0;

    for (auto &poly_line: line_strings) {
        total_segs += poly_line->getCoordinates()->size() - 1;
        base_vec.push_back(std::move(poly_line));
        base_ptrs.push_back(base_vec.back().get());
//        auto coords = poly_line->getCoordinates();
//
//        for (size_t i = 0; i < coords->size() - 1; i++) {
//            auto x1 = coords->getAt(i).x;
//            auto y1 = coords->getAt(i).y;
//            auto x2 = coords->getAt(i + 1).x;
//            auto y2 = coords->getAt(i + 1).y;
//            geos::geom::LineSegment l(x1, y1, x2, y2);
//
//            base_vec.push_back(l.toGeometry(*global_factory));
//            base_ptrs.push_back(base_vec.back().get());
//        }
    }
    Stopwatch sw;
    sw.start();
    index.glin_bulk_load(base_ptrs, error_bound, "z", -180.0, -180.0, 0.0000005, 0.0000005, pieces);
    sw.stop();

    printf("Bulk load %f\n", sw.ms());

    size_t total_xsects = 0;
    sw.start();

#pragma omp parallel
    {
#pragma omp  for reduction(+:total_xsects)
        for (auto &q: queries) {
            auto coords = q->getCoordinates();
            std::vector<geos::geom::Geometry *> find_result;
            double query_min_x = std::numeric_limits<double>::max();
            double query_max_x = std::numeric_limits<double>::lowest();
            double query_min_y = std::numeric_limits<double>::max();
            double query_max_y = std::numeric_limits<double>::lowest();

            for (size_t i = 0; i < coords->size() - 1; i++) {
                auto x1 = coords->getAt(i).x;
                auto y1 = coords->getAt(i).y;
                auto x2 = coords->getAt(i + 1).x;
                auto y2 = coords->getAt(i + 1).y;
                auto min_x = std::min(x1, x2);
                auto max_x = std::max(x1, x2);
                auto min_y = std::min(y1, y2);
                auto max_y = std::max(y1, y2);

                query_min_x = std::min(query_min_x, min_x);
                query_max_x = std::max(query_max_x, max_x);
                query_min_y = std::min(query_min_y, min_y);
                query_max_y = std::max(query_max_y, max_y);
            }
            auto *query = CreateQuery(query_min_x, query_min_y, query_max_x, query_max_y, global_factory);
            int count_filter = 0;
            find_result.clear();
            index.glin_find(query, "z", -180.0, -180.0, 0.0000005, 0.0000005, pieces, find_result, count_filter);

            for (size_t i = 0; i < coords->size() - 1; i++) {
                auto x1 = coords->getAt(i).x;
                auto y1 = coords->getAt(i).y;
                auto x2 = coords->getAt(i + 1).x;
                auto y2 = coords->getAt(i + 1).y;
                geos::geom::LineSegment l(x1, y1, x2, y2);
                for (auto *g: find_result) {
                    total_xsects += CalculateIntersections(g, l);
                }
            }
        }
    }
    sw.stop();

    printf("Total xsect: %zu\n", total_xsects.load());
    printf("Search time %f\n", sw.ms());
}
