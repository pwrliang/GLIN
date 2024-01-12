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

int CalculateIntersections(geos::geom::Geometry *g,
                           const geos::geom::LineSegment &l2) {
    auto seq = g->getCoordinates();
    int n_xsects = 0;

    for (size_t i = 0; i < seq->size() - 1; i++) {

        geos::geom::LineSegment l1(seq->getAt(i), seq->getAt(i + 1));
        auto c = l1.intersection(l2);
        if (!c.isNull()) {
            n_xsects++;
        }
    }
    return n_xsects;
}


#endif //GLIN_READ_WKT_H
