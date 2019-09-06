/* This file is part of the EpiCell library.
*
* The most recent release of EpiCell can be downloaded at
* <http://epicells.unige.ch/>
*
* EpiCell is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* EpiCell is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with EpiCell.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DATA_FUNCTIONALS_H
#define DATA_FUNCTIONALS_H

#include "mathematics/constants.h"
#include "core/cellCluster.h"
#include "environment/elasticEnvironment.h"
#include "environment/stickyEnvironment.h"
#include "environment/gradientEnvironment.h"
#include <cmath>

namespace epc {

namespace dataFunctionals2D {

// ============== cellCluster functionals ================== //

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getId = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->getId());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> computeArea = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->computeArea());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> computeApicoBasalOverIntercellularAspectRatio = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->computeApicoBasalOverIntercellularAspectRatio());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> computeAge = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back((double)c->getBirthDate());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> computeDegree = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->computeDegree());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> isOnMitosis = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->isOnMitosis());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getContainingVertex = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->getContainingVertex());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getType = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->getType().getId());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getRankToSource = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->getRankToSource());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getReferenceSourceCell = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        auto refCell = c->getReferenceSourceCell();
        if (refCell){
            vals.push_back(c->getReferenceSourceCell()->getId());
        } else {
            vals.push_back(c->getId());
        }
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getDistanceToSource = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->computeDistanceToSource());
    }
    return vals;
};


std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> computePositionAlongMajorAxis = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    std::vector<std::shared_ptr<Cell<2>>> DPCells;
    for (auto c: cells){
        if (c->getRankToSource()==0){
            DPCells.push_back(c);
        }
    }
    Array<double,2> centerDP = dataAnalysis::computeCentroid(DPCells);
    // Compute eignevalues and eigenvectors.
    std::vector<std::tuple<Array<double, 2>, double>> eigenVectorsValues = dataAnalysis::computePrincipleAxes(DPCells);
    Array<double,2> maxEigVect = std::get<0>(eigenVectorsValues[0]);
    for (auto c:cells){
        vals.push_back(c->computeProjectionOnAxis(centerDP, maxEigVect));
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> computePositionAlongMinorAxis = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    std::vector<std::shared_ptr<Cell<2>>> DPCells;
    for (auto c: cells){
        if (c->getRankToSource()==0){
            DPCells.push_back(c);
        }
    }
    Array<double,2> centerDP = dataAnalysis::computeCentroid(DPCells);
    // Compute eignevalues and eigenvectors.
    std::vector<std::tuple<Array<double, 2>, double>> eigenVectorsValues = dataAnalysis::computePrincipleAxes(DPCells);
    Array<double,2> minEigVect = std::get<0>(eigenVectorsValues[1]);
    for (auto c:cells){
        vals.push_back(c->computeProjectionOnAxis(centerDP, minEigVect));
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getProlifSignal = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->getSignal());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getMolE = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(c->getMolE());
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getSqrProlifSignal = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(pow(c->getSignal(),2));
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> getCubicProlifSignal = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells){
    std::vector<double> vals;
    for (auto c: cells) {
        vals.push_back(pow(c->getSignal(),3));
    }
    return vals;
};

std::function<std::vector<Array<double,2>> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> computeForce = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells) {
    std::vector<Array<double,2>> vals;
    for (auto c: cells) {
        const auto verts = vertexHelpers<2>::sortVerticesCounterClockwise(cellHelpers<2>::getVertices(*c));
        for (auto v: verts) {
            vals.push_back(v->getForce());
        }
    }
    return vals;
};

std::function<std::vector<double> (const std::vector<std::shared_ptr<Cell<2>>> &cells)> isInCell = [] (const std::vector<std::shared_ptr<Cell<2>>> &cells) {
    std::vector<double> vals;
    for (auto c: cells) {
        const auto verts = vertexHelpers<2>::sortVerticesCounterClockwise(cellHelpers<2>::getVertices(*c));
        for (auto v: verts) {
            vals.push_back(v->getInCell());
        }
    }
    return vals;
};

// ============== force functionals ================== //

auto constForce = [](const Vertex<2> &v, const Array<double,2> &f) -> Array<double,2> { return f;};

// ============== force functionals ================== //

auto radialForce = [](const Vertex<2> &v, const Array<double,2> &center, double amplitude) -> Array<double,2>
{
    Array<double,2> r = v.getPosition()-center;
    double distance = r.norm();
    r /= distance;

    return amplitude * r;
};

// ============== force functionals ================== //

auto gradientForce = [](const Vertex<2> &v, const GradientEnvironment<2> &gradientEnvironment) -> Array<double,2>
{
    Array<double,2> r = gradientEnvironment.computeForceOnVertex(v);
    // std::cout << "vertex = " << &v << '\n';
    // std::cout << "force = " << r << "\n" << '\n';

    return r;
};

auto elasticForce = [](const Vertex<2> &v, const ElasticEnvironment<2> &elasticEnvironment) -> Array<double,2>
{
    Array<double,2> r = elasticEnvironment.computeForceOnVertex(v);

    return r;
};


auto stickyForce = [](const Vertex<2> &v, const StickyEnvironment<2> &stickyEnvironment) -> Array<double,2>
{
    Array<double,2> r = stickyEnvironment.computeForceOnVertex(v);

    return r;
};

auto adherenceForce = [](const Edge<2> &e, const StickyEnvironment<2> &stickyEnvironment) -> Array<double,2>
{
    // IMPORTANT! This returns the force applied on e->getVertexOne() = (-1)* force on e->getVertexTwo().
    return stickyEnvironment.computeForceAlongEdge(e);
};

} // end dataFunctionals2D


}  // end namespace epc

#endif
