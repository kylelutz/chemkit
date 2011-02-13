#version 120

/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

varying vec3 vertex;
varying vec3 normal;

void main()
{
    vec3 L = normalize(gl_LightSource[0].position.xyz - vertex);
    vec3 E = normalize(-vertex);
    vec3 R = normalize(-reflect(L, normal));

    vec4 ambient = gl_FrontLightProduct[0].ambient;
    vec4 diffuse = gl_FrontLightProduct[0].diffuse * max(dot(normal, L), 0.0);
    vec4 specular = gl_FrontLightProduct[0].specular * pow(max(dot(R, E), 0.0), gl_FrontMaterial.shininess);

    gl_FragColor = gl_FrontLightModelProduct.sceneColor + ambient + diffuse + specular;
}
