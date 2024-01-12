#ifndef GLIN_READ_WKT_H
#define GLIN_READ_WKT_H

#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <fstream>
#include <geos/geom/util/PolygonExtracter.h>

/* For WKT read/write */
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>

std::vector<std::unique_ptr<geos::geom::Geometry>> ReadWKT(const std::string &path) {
    /* Geometry/GeometryFactory */
    using namespace geos::geom;

/* WKTReader/WKTWriter */
    using namespace geos::io;
    /* New factory with default (float) precision model */
    GeometryFactory::Ptr factory = GeometryFactory::create();

    /*
    * Reader requires a factory to bind the geometry to
    * for shared resources like the PrecisionModel
    */
    WKTReader reader(*factory);

    std::ifstream ifs(path);
    std::string line;
    std::vector<std::unique_ptr<geos::geom::Geometry>> geos;

    while (std::getline(ifs, line)) {
        if (!line.empty()) {
            std::unique_ptr<Geometry> poly = reader.read(line);
            geos.push_back(std::move(poly));
        }
    }
    return geos;
}

std::vector<geos::geom::LineSegment> GeometryToLineSegs(geos::geom::Geometry *g) {
    std::vector<geos::geom::LineSegment> segs;
    auto seq = g->getCoordinates();

    for (size_t i = 0; i < seq->size() - 1; i++) {

        geos::geom::LineSegment l(seq->getAt(i), seq->getAt(i + 1));
        segs.push_back(l);
    }

    return segs;
}

//std::vector<std::unique_ptr<geos::geom::LineString>>
//PolygonsToLineSegs(std::vector<std::unique_ptr<geos::geom::Geometry>> &polygons) {
//    std::vector<std::unique_ptr<geos::geom::LineString>> segs;
//
//    for (auto &geo: polygons) {
//        std::vector<const geos::geom::Polygon *> polys;
//        geos::geom::util::PolygonExtracter::getPolygons(*geo, polys);
//
//        for (auto *p: polys) {
//            auto seq = p->getExteriorRing()->getCoordinates();
//
//            for (size_t i = 0; i < seq->size() - 1; i++) {
//                geos::geom::LineSegment l(seq->getAt(i), seq->getAt(i + 1));
//                segs.push_back(l.toGeometry(*factory));
//            }
//        }
//    }
//    return segs;
//}

std::vector<geos::geom::LineSegment>
PolygonsToLineSegs(std::vector<std::unique_ptr<geos::geom::Geometry>> &polygons) {
    std::vector<geos::geom::LineSegment> segs;

    for (auto &geo: polygons) {
        std::vector<const geos::geom::Polygon *> polys;
        geos::geom::util::PolygonExtracter::getPolygons(*geo, polys);

        for (auto *p: polys) {
            auto seq = p->getExteriorRing()->getCoordinates();

            for (size_t i = 0; i < seq->size() - 1; i++) {
                geos::geom::LineSegment l(seq->getAt(i), seq->getAt(i + 1));
                segs.push_back(l);
            }
        }
    }
    return segs;
}

struct edge_params {
    double a, b, c;

    void from(const geos::geom::LineSegment &l) {
        auto x1 = l.p0.x;
        auto y1 = l.p0.y;
        auto x2 = l.p1.x;
        auto y2 = l.p1.y;

        a = y1 - y2;
        b = x2 - x1;
        c = -x1 * a - y1 * b;

        assert(a != 0 || b != 0);

        if (b < 0) {
            a = -a;
            b = -b;
            c = -c;
        }
    }
};

int CalculateIntersections(const geos::geom::LineSegment &l1,
                           const geos::geom::LineSegment &l2) {
    edge_params e1, e2;
    auto e1_p1 = l1.p0, e1_p2 = l1.p1;
    auto e2_p1 = l2.p0, e2_p2 = l2.p1;
    e1.from(l1);
    e2.from(l2);

#define SUBEDGE(p, e) \
  (p.x * e.a + p.y * e.b + e.c)  // ax+by+c
    // N.B., e.b >= 0, we ensure it when calculate edge eqns
    auto e2_p1_agst_e1 = SUBEDGE(e2_p1, e1);
    auto e2_p2_agst_e1 = SUBEDGE(e2_p2, e1);
    auto e1_p1_agst_e2 = SUBEDGE(e1_p1, e2);
    auto e1_p2_agst_e2 = SUBEDGE(e1_p2, e2);
#undef SUBEDGE
    // e1_p1 is on e2
    if (e1_p1_agst_e2 == 0) {
        e1_p1_agst_e2 = -e2.a;
    }
    if (e1_p1_agst_e2 == 0) {  // a = 0, e2 is parallel to x-axis
        e1_p1_agst_e2 = -e2.b;
    }
    if (e1_p1_agst_e2 == 0) {  // b = 0, then c must be 0
        return false;            // zero length edge
    }

    if (e1_p2_agst_e2 == 0) {
        e1_p2_agst_e2 = -e2.a;
    }
    if (e1_p2_agst_e2 == 0) {
        e1_p2_agst_e2 = -e2.b;
    }
    if (e1_p2_agst_e2 == 0) {
        return false;
    }

    // p1 and p2 of edge1 is on the same side of edge2, they will not
    // intersect
    if ((e1_p1_agst_e2 > 0 && e1_p2_agst_e2 > 0) ||
        (e1_p1_agst_e2 < 0 && e1_p2_agst_e2 < 0)) {
        return false;
    }

    // e2_p1 is on e1
    if (e2_p1_agst_e1 == 0) {
        e2_p1_agst_e1 = e1.a;
    }
    if (e2_p1_agst_e1 == 0) {
        e2_p1_agst_e1 = e1.b;
    }
    if (e2_p1_agst_e1 == 0) {
        return false;
    }
    if (e2_p2_agst_e1 == 0) {
        e2_p2_agst_e1 = e1.a;
    }
    if (e2_p2_agst_e1 == 0) {
        e2_p2_agst_e1 = e1.b;
    }
    if (e2_p2_agst_e1 == 0) {
        return false;
    }
    if ((e2_p1_agst_e1 > 0 && e2_p2_agst_e1 > 0) ||
        (e2_p1_agst_e1 < 0 && e2_p2_agst_e1 < 0)) {
        return false;
    }

    /*
     * Check if both edges are the same.  If so, they shouldn't be
     * intersecting.
     */
    if ((e1_p1 == e2_p1 && e1_p2 == e2_p2) ||
        (e1_p1 == e2_p2 && e1_p2 == e2_p1)) {
        return false;
    }

    return true;
}

int CalculateIntersections(geos::geom::Geometry *g,
                           const geos::geom::LineSegment &l2) {
    auto seq = g->getCoordinates();
    int n_xsects = 0;

    for (size_t i = 0; i < seq->size() - 1; i++) {
        geos::geom::LineSegment l1(seq->getAt(i), seq->getAt(i + 1));
        if (CalculateIntersections(l1, l2)) {
            n_xsects++;
        }
    }
    return n_xsects;
}


#endif //GLIN_READ_WKT_H
