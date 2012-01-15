/******************************************************************************
**
** Copyright (C) 2011 Wang Wenlin <sopl.wang@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#ifndef CHEMKIT_GRAPHICSSOLVENTSURFACEITEM_H
#define CHEMKIT_GRAPHICSSOLVENTSURFACEITEM_H

#include "graphics.h"

#include "graphicsitem.h"

namespace chemkit {

class Molecule;
class GraphicsSolventSurfaceItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsSolventSurfaceItem : public GraphicsItem
{
public:
    // enumerations
    enum SurfaceQuality {
        SurfaceQualityMoreMiserable = -4,
        SurfaceQualityMiserable = -3,
        SurfaceQualityMorePoor = -2,
        SurfaceQualityPoor = -1,
        SurfaceQualityNormal = 0,
        SurfaceQualityGood = 1,
        SurfaceQualityNearPerfect = 2,
        SurfaceQualityPerfect = 3,
        SurfaceQualityImpractical = 4
    };

    enum SurfaceType {
        SurfaceTypeSolid = 0,
        SurfaceTypeDots = 1,
        SurfaceTypeTriangles = 2,
        SurfaceTypeType3 = 3,
        SurfaceTypeType4 = 4,
        SurfaceTypeType5 = 5,
        SurfaceTypeType6 = 6
    };

    enum SurfaceSolventType {
        SurfaceSolventTypeExcluded = 0,
        SurfaceSolventTypeAccessible = 1
    };

    // construction and destruction
    GraphicsSolventSurfaceItem(const Molecule *molecule = 0, SurfaceSolventType type = SurfaceSolventTypeExcluded);
    ~GraphicsSolventSurfaceItem();

    // properties
    void setMolecule(const Molecule *molecule);
    const Molecule* molecule() const;
    void setQuality(SurfaceQuality quality);
    SurfaceQuality quality() const;
    void setSurfaceType(SurfaceType type);
    SurfaceType surfaceType() const;
    void setSurfaceSolventType(SurfaceSolventType type);
    SurfaceSolventType surfaceSolventType() const;
    void setProbeRadius(Real radius);
    Real probeRadius() const;
    void setColor(const QColor &color);
    QColor color() const;

    // drawing
    virtual void paint(GraphicsPainter *painter);

private:
    // internal methods
    void setCalculated(bool calculated);
    Real maxVdwRadius();

private:
    GraphicsSolventSurfaceItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSSOLVENTSURFACEITEM_H
