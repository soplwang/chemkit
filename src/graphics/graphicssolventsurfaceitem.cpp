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

#include "graphicssolventsurfaceitem.h"

#include <chemkit/geometry.h>
#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <algorithm>

#include "graphicspainter.h"
#include "graphicsvertexbuffer.h"

extern "C" {
#include "../3rdparty/mskit/MSKContext.h"
#include "../3rdparty/mskit/MemoryDebug.h"
#include "../3rdparty/mskit/SurfaceJob.h"
}

namespace chemkit {

namespace {

class __mskit_context_helper {
public:
    __mskit_context_helper() : ctx(MSKContextNew()) {}
    ~__mskit_context_helper() { MSKContextFree(ctx); }

public:
    MSKContext *ctx;
};

GraphicsVertexBuffer* calculateSurface(const std::vector<Point3>& points,
                                       const std::vector<Real>& radii,
                                       Real max_vdw, Real probe_radius,
                                       int surface_quality, int surface_type,
                                       int surface_solvent)
{
    static __mskit_context_helper _ctx_holder;

    float *coord = VLAlloc(float, points.size() * 3);
    SurfaceJobAtomInfo *atom_info = VLACalloc(SurfaceJobAtomInfo, points.size());

    if (coord && atom_info) {
        float *cp = coord;
        SurfaceJobAtomInfo *ap = atom_info;

        for (std::vector<Point3>::const_iterator i = points.begin(); i < points.end(); i++) {
            *cp++ = i->x();
            *cp++ = i->y();
            *cp++ = i->z();
        }
        for (std::vector<Real>::const_iterator i = radii.begin(); i < radii.end(); i++) {
            (ap++)->vdw = static_cast<float>(*i);
        }

        SurfaceJob *job = SurfaceJobNew(_ctx_holder.ctx, coord, atom_info,
                                        max_vdw, probe_radius,
                                        surface_quality, surface_type, surface_solvent,
                                        10, 0, 7.0F,
                                        -3.0F, 0.2F, 2.0F);

        if (job && SurfaceJobRun(_ctx_holder.ctx, job)) {
            QVector<Point3f> verticies(job->N);
            QVector<Vector3f> normals(job->N);
            QVector<unsigned short> indicies(job->NT);

            fprintf(stderr, "job->N = %d, job->NT = %d\n", job->N, job->NT);

            for (float *vp = job->V, *np = job->VN, *e = (job->V + job->N*3); vp < e; vp+=3, np+=3) {
                verticies.push_back(Point3f(vp[0], vp[1], vp[2]));
                normals.push_back(Point3f(np[0], np[1], np[2]));
                fprintf(stderr, "v = %f,%f,%f  n = %f,%f,%f\n", vp[0], vp[1], vp[2], np[0], np[1], np[2]);
            }
            for (int *tp = job->T, *e = (job->T + job->NT); tp < e; tp++) {
                indicies.push_back(static_cast<unsigned short>(*tp));
                fprintf(stderr, "%d ", *tp);
            }

            // create vertex buffer
            GraphicsVertexBuffer *buffer = new GraphicsVertexBuffer;

            buffer->setVerticies(verticies);
            buffer->setNormals(normals);
            //buffer->setIndicies(indicies);

            SurfaceJobFree(_ctx_holder.ctx, job);

            return buffer;
        }

        if (job)
            SurfaceJobFree(_ctx_holder.ctx, job);
    }

    FreeP(atom_info);
    FreeP(coord);

    return 0;
}

}

// === GraphicsSolventSurfaceItemPrivate ======================================= //
class GraphicsSolventSurfaceItemPrivate
{
public:
    const Molecule *molecule;
    GraphicsSolventSurfaceItem::SurfaceQuality quality;
    GraphicsSolventSurfaceItem::SurfaceType surfaceType;
    GraphicsSolventSurfaceItem::SurfaceSolventType surfaceSolventType;
    Real probeRadius;
    QColor color;
    std::vector<Point3> points;
    std::vector<Real> radii;
    Real maxVdwRadius;
    bool maxVdwCalculated;
    GraphicsVertexBuffer *buffer;
};

// === GraphicsSolventSurfaceItem ============================================== //
/// \class GraphicsSolventSurfaceItem graphicssolventsurfaceitem.h chemkit/graphicssolventsurfaceitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsSolventSurfaceItem class visually displays a
///        solvent surface.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new solvent surface item to display of \p molecule.
GraphicsSolventSurfaceItem::GraphicsSolventSurfaceItem(const Molecule *molecule, SurfaceSolventType type)
    : GraphicsItem(),
      d(new GraphicsSolventSurfaceItemPrivate)
{
    d->molecule = molecule;
    d->quality = SurfaceQualityNormal;
    d->surfaceType = SurfaceTypeSolid;
    d->surfaceSolventType = type;
    d->probeRadius = 1.4;

    if (molecule) {
        d->points.reserve(molecule->size());
        d->radii.reserve(molecule->size());

        for(size_t i = 0; i < molecule->size(); i++){
            const Atom *atom = molecule->atom(i);

            d->points.push_back(atom->position());
            d->radii.push_back(atom->vanDerWaalsRadius());
        }
    }

    d->maxVdwCalculated = false;
    d->buffer = 0;
}

/// Destroys the solvent surface object.
GraphicsSolventSurfaceItem::~GraphicsSolventSurfaceItem()
{
    delete d->buffer;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule for the surface.
void GraphicsSolventSurfaceItem::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;

    // update atom positions and radii
    if(molecule){
        d->points.resize(molecule->size());
        d->radii.resize(molecule->size());

        for(size_t i = 0; i < molecule->size(); i++){
            const Atom *atom = molecule->atom(i);

            d->points[i] = atom->position();
            d->radii[i] = atom->vanDerWaalsRadius();
        }
    }

    setCalculated(false);
}

/// Returns the molecule for the surface.
const Molecule* GraphicsSolventSurfaceItem::molecule() const
{
    return d->molecule;
}

/// Sets the surface quality to \p quality.
void GraphicsSolventSurfaceItem::setQuality(SurfaceQuality quality)
{
    d->quality = quality;
    setCalculated(false);
}

/// Returns the surface quality.
GraphicsSolventSurfaceItem::SurfaceQuality GraphicsSolventSurfaceItem::quality() const
{
    return d->quality;
}

/// Sets the surface type to \p type.
void GraphicsSolventSurfaceItem::setSurfaceType(SurfaceType type)
{
    d->surfaceType = type;
    setCalculated(false);
}

/// Returns the surface type.
GraphicsSolventSurfaceItem::SurfaceType GraphicsSolventSurfaceItem::surfaceType() const
{
    return d->surfaceType;
}

/// Sets the surface solvent type to \p type.
void GraphicsSolventSurfaceItem::setSurfaceSolventType(SurfaceSolventType type)
{
    d->surfaceSolventType = type;
    setCalculated(false);
}

/// Returns the surface solvent type.
GraphicsSolventSurfaceItem::SurfaceSolventType GraphicsSolventSurfaceItem::surfaceSolventType() const
{
    return d->surfaceSolventType;
}

/// Sets the probe radius to \p radius.
void GraphicsSolventSurfaceItem::setProbeRadius(Real radius)
{
    d->probeRadius = radius;
    setCalculated(false);
}

/// Returns the probe radius.
///
/// The default probe radius is 1.4 Angstroms which approximates
/// the radius of a water molecule.
Real GraphicsSolventSurfaceItem::probeRadius() const
{
    return d->probeRadius;
}

/// Sets the color for the solvent surface.
void GraphicsSolventSurfaceItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the color for the solvent surface.
QColor GraphicsSolventSurfaceItem::color() const
{
    return d->color;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsSolventSurfaceItem::paint(GraphicsPainter *painter)
{
    if(!d->molecule){
        return;
    }

    if(!d->buffer){
        d->buffer = calculateSurface(d->points, d->radii,
                                     maxVdwRadius(), d->probeRadius,
                                     d->quality, d->surfaceType, d->surfaceSolventType);
        if(!d->buffer) {
            return;
        }
    }

    QColor color = d->color;
    color.setAlphaF(opacity());

    painter->setColor(color);
    painter->draw(d->buffer);
}

// --- Internal Methods ---------------------------------------------------- //
void GraphicsSolventSurfaceItem::setCalculated(bool calculated)
{
    if (!calculated) {
        delete d->buffer;
        d->buffer = 0;
        d->maxVdwCalculated = false;
    }
}

Real GraphicsSolventSurfaceItem::maxVdwRadius()
{
    if (!d->maxVdwCalculated) {
        if (d->radii.empty()) {
            d->maxVdwRadius = 0;
        } else {
            d->maxVdwRadius = *std::max_element(d->radii.begin(), d->radii.end());
        }
        d->maxVdwCalculated = true;
    }
    return d->maxVdwRadius;
}

} // end chemkit namespace
