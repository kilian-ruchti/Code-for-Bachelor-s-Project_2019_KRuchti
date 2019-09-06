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

#ifndef GRADIENTENVIRONMENT_HH
#define CELLTYPE_HH

#include "gradientEnvironment.h"
#include <math.h>

namespace epc {

template<int d>
GradientEnvironment<d>::GradientEnvironment(bool activatedFunction_, double xmin_, double xmax_, double ymin_, double ymax_) : activatedFunction(activatedFunction_), xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_){}


template<int d>
Array<double, d> GradientEnvironment<d>::computeForceOnVertex(const Vertex<d> &v) const{

  if (activatedFunction){
		// std::cout << "activatedFunction at TRUE !!" << std::endl;
    return forceApplied(v.getPosition());
	// 	return elasticCoef*computeDistanceFromBorders(v.getPosition())*computeForceDirection(v.getPosition());
	}
	return Array<double, d>(double(0.0));
}

template<int d>
Array<double, d> GradientEnvironment<d>::computeForceDirection(Array<double, d> position) const{
	Array<double, d> ref(double(0.0));
	ref[0] = 1;
	// std::cout << "GradientEnvironment::computeForceDirection(...) = " << ref << std::endl;
	return ref;
}

// template<int d>
// Array<double, d> GradientEnvironment<d>::forceApplied(Array<double, d> pos) const{
// // void GradientEnvironment<d>::forceApplied(Array<double, d> position) const{
//
//   // Array<double,2> force(1e-7,0.0);
//   Array<double,2> force(2e-9, 0.0);
//   double xmin = getXmin();
//   double xmax = getXmax();
//
//   double a = xmin+(xmax-xmin)/3;
//   double b = a+(xmax-xmin)/3;
//   double h = (b-a)/10;
//
//   if (pos[0] < a) {
//     force[0] = 0;
//   } else if (pos[0] < a+h) {
//     force[0] = force[0]*(10*pos[0]);
//   } else  if (pos[0] < a+(3*h)) {
//     force[0] = force[0]*(-2*pos[0]+12);
//   } else if (pos[0] < a+(6*h)) {
//     force[0] = force[0]*(-((4*pos[0])/3) + 10);
//   } else if (pos[0] < b) {
//     force[0] = force[0]*(-(pos[0]/2)+5);
//   } else {
//     force[0] = 0;
//   }
//
// 	return force;
// }


// template<int d>
// Array<double, d> GradientEnvironment<d>::forceApplied(Array<double, d> pos) const{
// // void GradientEnvironment<d>::forceApplied(Array<double, d> position) const{
//
//   Array<double,2> force(1e-8, 0.0);
//   double xmin = getXmin();
//   double xmax = getXmax();
//
//   double h = (xmax-xmin)/10;
//
//   if (pos[0] < xmin+h) {
//     force[0] = force[0]*(-2/h * pos[0] + (2*xmin/h + 4));
//   } else if (pos[0] < xmin+4*h) {
//     force[0] = force[0]*2;
//   } else  if (pos[0] < xmin+5*h) {
//     force[0] = force[0]*(-2/h * pos[0] + (2*xmin/h + 10));
//   } else if (pos[0] < xmin+6*h) {
//     force[0] = force[0]*(8/h * pos[0] - (8*xmin/h + 40));
//   } else if (pos[0] < xmin+(25/4)*h) {
//     force[0] = force[0]*(-28/h * pos[0] + (28*xmin/h + 176));
//   } else if (pos[0] < xmin+(26/4)*h) {
//     force[0] = force[0]*(12/h * pos[0] - (12*xmin/h + 74));
//   } else if (pos[0] < xmin+7*h) {
//     force[0] = force[0]*(-4/h * pos[0] + (4*xmin/h + 30));
//   } else if (pos[0] > xmin+7*h) {
//     force[0] = force[0]*2;
//   }
// 	return force;
// }

template<int d>
Array<double, d> GradientEnvironment<d>::forceApplied(Array<double, d> pos) const{
// void GradientEnvironment<d>::forceApplied(Array<double, d> position) const{

  // Array<double,2> force(1e-7,0.0);
  Array<double,2> force(1e-8, 0.0);
  double xmin = getXmin();
  double xmax = getXmax();

  double h = (xmax-xmin)/10;

  if (pos[0] < xmin+h) {
    force[0] = force[0]*(-2/h * pos[0] + (2*xmin/h + 4));
  } else if (pos[0] < xmin+4*h) {
    force[0] = force[0]*2;
  } else  if (pos[0] < xmin+5*h) {
    force[0] = force[0]*(-2/h * pos[0] + (2*xmin/h + 10));
  } else if (pos[0] < xmin+6*h) {
    force[0] = force[0]*(8/h * pos[0] - (8*xmin/h + 40));
  } else if (pos[0] < xmin+7*h) {
    force[0] = 0;
  } else if (pos[0] > xmin+7*h) {
    force[0] = force[0]*2;
  }
	return force;
}


template<int d>
double const& GradientEnvironment<d>::getXmin() const {
  return xmin;
};

template<int d>
double const& GradientEnvironment<d>::getXmax() const {
  return xmax;
};

template<int d>
double const& GradientEnvironment<d>::getYmin() const {
  return ymin;
};

template<int d>
double const& GradientEnvironment<d>::getYmax() const {
  return ymax;
};

}  // end namespace epc

#endif
