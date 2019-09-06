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

// #ifndef GRADIENTENVIRONMENT_H
// #define GRADIENTENVIRONMENT_H
//
// #include <core/vertex.h>
// #include <mathematics/array.h>
//
// namespace epc {
//
// template<int d>
// class GradEnvironment {
//
// public:
// 	Array<double, d> computeForceOnVertex(const Vertex<d> &v) const;
// };
//
//
// template<int d>
// class GradientEnvironment : public GradEnvironment<d>{
//
// public:
// 	GradientEnvironment(const bool activatedFunction_, double xmin_, double xmax_, double ymin_, double ymax_);
// 	Array<double, d> computeForceOnVertex(const Vertex<d> &v) const;
// private:
// 	Array<double, d> computeForceDirection(Array<double, d> position) const;
//
// private:
// 	const bool activatedFunction;
// 	double xmin;
// 	double xmax;
// 	double ymin;
// 	double ymax;
// };
// }   // end namespace epc
//
// #endif


#ifndef GRADIENTENVIRONMENT_H
#define GRADIENTENVIRONMENT_H

#include <core/vertex.h>
#include <mathematics/array.h>

namespace epc {

template<int d>
class GradientEnvironment {

public:
	GradientEnvironment(bool activatedFunction_, double xmin_, double xmax_, double ymin_, double ymax_);
	Array<double, d> computeForceOnVertex(const Vertex<d> &v) const;
private:
	Array<double, d> computeForceDirection(Array<double, d> position) const;
	Array<double, d> forceApplied(Array<double, d> position) const;
	// void forceApplied(Array<double, d> position) const;
	double const& getXmin() const;
	double const& getXmax() const;
	double const& getYmin() const;
	double const& getYmax() const;

private:
	bool activatedFunction;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
};
}   // end namespace epc

#endif
