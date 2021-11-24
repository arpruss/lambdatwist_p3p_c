/*
P3P
Copyright (C) 2021 Mikael Persson, Klas Nordberg and Alexander Pruss

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _P3P_H
#define _P3P_H

typedef double Vector2[2];
typedef double Vector3[3];
typedef double Matrix33[3][3];

int p3p_lambdatwist(Vector2 cam1, Vector2 cam2, Vector2 cam3, Vector3 x1, Vector3 x2, Vector3 x3,
        Matrix33* Rs, Vector3* Ts, int refinement_iterations);
        
void toCam(Vector2 out, Vector3 in,Matrix33 R,Vector3 T) {
    Vector3 temp;
    apply3(temp, R, in);
    add3(temp, T, temp);
    out[0] = temp[0]/temp[2];
    out[1] = temp[1]/temp[2];
}

#endif // _P3P_H       