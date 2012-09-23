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

#include "graphicspymolsurfaceitem.h"

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <chemkit/atom.h>
#include <chemkit/element.h>
#include <chemkit/foreach.h>
#include <chemkit/geometry.h>
#include <chemkit/molecule.h>

#include "graphicspainter.h"
#include "graphicsmaterial.h"
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
                                       const std::vector<Element>& elements,
                                       Real max_vdw, Real probe_radius,
                                       int surface_quality, int surface_type, int surface_solvent,
                                       const AtomColorMap& colorMap, float opacity)
{
    static __mskit_context_helper _ctx_holder;

    MSKContextClean(_ctx_holder.ctx);

    float *coord = VLAlloc(float, points.size() * 3);
    SurfaceJobAtomInfo *atom_info = VLACalloc(SurfaceJobAtomInfo, points.size());

    if (coord && atom_info) {
        float *cp = coord;
        SurfaceJobAtomInfo *ap = atom_info;

        for (std::vector<Point3>::const_iterator i = points.begin(); i < points.end(); ++i) {
            *cp++ = static_cast<float>(i->x());
            *cp++ = static_cast<float>(i->y());
            *cp++ = static_cast<float>(i->z());
        }
        for (std::vector<Real>::const_iterator i = radii.begin(); i < radii.end(); ++i) {
            (ap++)->vdw = static_cast<float>(*i);
        }

        SurfaceJob *job = SurfaceJobNew(_ctx_holder.ctx, coord, atom_info,
                                        max_vdw, probe_radius,
                                        surface_quality, surface_type, surface_solvent,
                                        10, 0, 7.0F, -3.0F,
                                        0.2F, 2.0F);

        if (job && SurfaceJobRun(_ctx_holder.ctx, job)) {
            QVector<Point3f> vertices;
            QVector<Vector3f> normals;
            QVector<unsigned short> indices;

            vertices.reserve(job->N);
            normals.reserve(job->N);
            for (float *vp = job->V, *np = job->VN, *e = (job->V + job->N*3); vp < e; vp+=3, np+=3) {
                vertices.push_back(Point3f(vp[0], vp[1], vp[2]));
                normals.push_back(Point3f(np[0], np[1], np[2]));
            }

            if (surface_type != 1) {
                indices.reserve(job->NT*3);

                for (int *tp = job->T, *e = (job->T + job->NT*3); tp < e; ++tp) {
                    indices.push_back(static_cast<unsigned short>(*tp));
                }
            }

            // create vertex buffer
            GraphicsVertexBuffer *buffer = new GraphicsVertexBuffer;

            buffer->setVertices(vertices);
            buffer->setNormals(normals);
            buffer->setIndices(indices);

            // apply colors
            if(!elements.empty() && SurfaceJobOwnership(_ctx_holder.ctx, job)) {
                QVector<QColor> colors;
                colors.reserve(job->N);

                for (int *op = job->VO, *e = (job->VO + job->N); op < e; ++op) {
                    QColor color = colorMap.color(elements[*op]);
                    color.setAlphaF(opacity);
                    colors.push_back(color);
                }

                buffer->setColors(colors);
            }

            SurfaceJobFree(_ctx_holder.ctx, job);

            return buffer;
        }

        if (job)
            SurfaceJobFree(_ctx_holder.ctx, job);
    }

    VLAFreeP(atom_info);
    VLAFreeP(coord);

    return 0;
}

// Returns the electrostatic potential at the given position calculated
// from the partial charges and positions of the atoms in the molecule.
float electrostaticPotential(const Molecule *molecule, const Point3f &position)
{
    float esp = 0;

    const float pi = chemkit::constants::Pi;
    const float e0 = 1.0f;

    foreach(const Atom *atom, molecule->atoms()){
        float q = atom->partialCharge();
        float r = (position - atom->position().cast<float>()).norm();

        esp += (1.0f / (4.0f * pi * e0)) * (q / r);
    }

    return esp;
}

// Returns a color interpolated at the given value between the colors a
// (starting at value av) and b (starting at value bv).
QColor interpolate(const QColor &a, const QColor &b, float av, float bv, float value)
{
    float v = (value - av) / (bv - av);

    return QColor::fromRgbF(a.redF() + (b.redF() - a.redF()) * v,
                            a.greenF() + (b.greenF() - a.greenF()) * v,
                            a.blueF() + (b.blueF() - a.blueF()) * v);
}

// Returns the color associated with the electrostatic potential.
QColor electrostaticPotentialColor(float esp)
{
    QColor red = QColor::fromRgb(255, 0, 0);
    QColor orange = QColor::fromRgb(255, 127, 0);
    QColor yellow = QColor::fromRgb(255, 255, 0);
    QColor green = QColor::fromRgb(0, 255, 0);
    QColor blue = QColor::fromRgb(0, 0, 255);

    // color ranges are hard coded (for now)
    float redStart = -0.0075f;
    float orangeStart = -0.0035f;
    float yellowStart = 0.0f;
    float greenStart = 0.0015f;
    float blueStart = 0.0045f;

    if(esp < redStart){
        return red;
    }
    else if(esp < orangeStart){
        return interpolate(red, orange, redStart, orangeStart, esp);
    }
    else if(esp < yellowStart){
        return interpolate(orange, yellow, orangeStart, yellowStart, esp);
    }
    else if(esp < greenStart){
        return interpolate(yellow, green, yellowStart, greenStart, esp);
    }
    else if(esp < blueStart){
        return interpolate(green, blue, greenStart, blueStart, esp);
    }
    else{
        return blue;
    }
}

} // end anonymous namespace

// === GraphicsPymolSurfaceItemPrivate ======================================= //
class GraphicsPymolSurfaceItemPrivate
{
public:
    const Molecule *molecule;
    GraphicsPymolSurfaceItem::SurfaceQuality quality;
    GraphicsPymolSurfaceItem::SurfaceType surfaceType;
    GraphicsPymolSurfaceItem::SolventType solventType;
    Real probeRadius;
    GraphicsPymolSurfaceItem::ColorMode colorMode;
    QColor color;
    boost::shared_ptr<AtomColorMap> colorMap;
    std::vector<Point3> points;
    std::vector<Real> radii;
    std::vector<Element> elements;
    Real maxVdwRadius;
    bool maxVdwCalculated;
    GraphicsVertexBuffer *buffer;
};

// === GraphicsPymolSurfaceItem ============================================ //
/// \class GraphicsPymolSurfaceItem graphicspymolsurfaceitem.h chemkit/graphicspymolsurfaceitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsPymolSurfaceItem class visually displays a Pymol style
///        solvent surface.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new solvent surface item to display of \p molecule.
GraphicsPymolSurfaceItem::GraphicsPymolSurfaceItem(const Molecule *molecule, SolventType solventType)
    : GraphicsItem(),
      d(new GraphicsPymolSurfaceItemPrivate)
{
    d->molecule = molecule;
    d->quality = SurfaceQualityNormal;
    d->surfaceType = SurfaceTypeSolid;
    d->solventType = solventType;
    d->probeRadius = 1.4;
    d->color = Qt::red;
    d->colorMode = AtomColor;
    d->colorMap = boost::make_shared<AtomColorMap>(AtomColorMap::DefaultColorScheme);

    if(molecule){
        d->points.reserve(molecule->size());
        d->radii.reserve(molecule->size());
        d->elements.reserve(molecule->size());

        foreach(const Atom *atom, molecule->atoms()){
            d->points.push_back(atom->position());
            d->radii.push_back(atom->vanDerWaalsRadius());
            d->elements.push_back(atom->element());
        }
    }

    d->maxVdwCalculated = false;
    d->buffer = 0;
}

/// Destroys the solvent surface object.
GraphicsPymolSurfaceItem::~GraphicsPymolSurfaceItem()
{
    delete d->buffer;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule for the surface.
void GraphicsPymolSurfaceItem::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;

    // update atom positions and radii
    if(molecule){
        d->points.resize(molecule->size());
        d->radii.resize(molecule->size());
        d->elements.resize(molecule->size());

        for(size_t i = 0; i < molecule->size(); i++){
            const Atom *atom = molecule->atom(i);

            d->points[i] = atom->position();
            d->radii[i] = atom->vanDerWaalsRadius();
            d->elements[i] = atom->element();
        }
    }

    setCalculated(false);
}

/// Returns the molecule for the surface.
const Molecule* GraphicsPymolSurfaceItem::molecule() const
{
    return d->molecule;
}

/// Sets the surface quality to \p quality.
void GraphicsPymolSurfaceItem::setQuality(SurfaceQuality quality)
{
    d->quality = quality;
    setCalculated(false);
}

/// Returns the surface quality.
GraphicsPymolSurfaceItem::SurfaceQuality GraphicsPymolSurfaceItem::quality() const
{
    return d->quality;
}

/// Sets the surface type to \p type.
void GraphicsPymolSurfaceItem::setSurfaceType(SurfaceType type)
{
    d->surfaceType = type;
    setCalculated(false);
}

/// Returns the surface type.
GraphicsPymolSurfaceItem::SurfaceType GraphicsPymolSurfaceItem::surfaceType() const
{
    return d->surfaceType;
}

/// Sets the surface solvent type to \p type.
void GraphicsPymolSurfaceItem::setSolventType(SolventType solventType)
{
    d->solventType = solventType;
    setCalculated(false);
}

/// Returns the surface solvent type.
GraphicsPymolSurfaceItem::SolventType GraphicsPymolSurfaceItem::solventType() const
{
    return d->solventType;
}

/// Sets the probe radius to \p radius.
void GraphicsPymolSurfaceItem::setProbeRadius(Real radius)
{
    d->probeRadius = radius;
    setCalculated(false);
}

/// Returns the probe radius.
///
/// The default probe radius is 1.4 Angstroms which approximates
/// the radius of a water molecule.
Real GraphicsPymolSurfaceItem::probeRadius() const
{
    return d->probeRadius;
}

/// Sets the color for the solvent surface.
void GraphicsPymolSurfaceItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the color for the solvent surface.
QColor GraphicsPymolSurfaceItem::color() const
{
    return d->color;
}

/// Sets the color mode for the solvent surface to \p mode.
void GraphicsPymolSurfaceItem::setColorMode(ColorMode mode)
{
    d->colorMode = mode;
    setCalculated(false);
}

/// Returns the color mode for the solvent surface.
GraphicsPymolSurfaceItem::ColorMode GraphicsPymolSurfaceItem::colorMode() const
{
    return d->colorMode;
}

/// Sets the color map for the solvent surface to \p colorMap.
void GraphicsPymolSurfaceItem::setColorMap(const boost::shared_ptr<AtomColorMap> &colorMap)
{
    d->colorMap = colorMap;
    setCalculated(false);
}

/// Returns the color map for the solvent surface.
boost::shared_ptr<AtomColorMap> GraphicsPymolSurfaceItem::colorMap() const
{
    return d->colorMap;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsPymolSurfaceItem::paint(GraphicsPainter *painter)
{
    if(!d->molecule){
        return;
    }

    if(!d->buffer){
        if(d->colorMode == SolidColor || d->colorMode == ElectrostaticPotential) {
            d->buffer = calculateSurface(d->points, d->radii, std::vector<Element>(),
                                         maxVdwRadius(), d->probeRadius,
                                         d->quality, d->surfaceType, d->solventType,
                                         *d->colorMap, 1.f);

            if(d->colorMode == ElectrostaticPotential){
                float opac = opacity();

                // function that returns electrostatic potential at a point in space
                boost::function<float (const Point3f&)> espFunction =
                    boost::bind(electrostaticPotential, d->molecule, _1);

                // calculate color for each vertex
                QVector<QColor> colors;
                colors.reserve(d->buffer->vertexCount());

                foreach(const Point3f &vertex, d->buffer->vertices()){
                    QColor color = electrostaticPotentialColor(espFunction(vertex));
                    color.setAlphaF(opac);
                    colors.append(color);
                }

                d->buffer->setColors(colors);
            }
        }
        else if(d->colorMode == AtomColor) {
            d->buffer = calculateSurface(d->points, d->radii, d->elements,
                                         maxVdwRadius(), d->probeRadius,
                                         d->quality, d->surfaceType, d->solventType,
                                         *d->colorMap, opacity());
        }

        if(!d->buffer) {
            return;
        }
    }

    if (d->colorMode == SolidColor) {
        QColor color = d->color;
        color.setAlphaF(opacity());
        painter->setColor(color);
    }

    painter->draw(d->buffer);
}

// --- Internal Methods ---------------------------------------------------- //
void GraphicsPymolSurfaceItem::itemChanged(ItemChange change)
{
    if(change == ItemOpacityChanged){
        if(isOpaque()){
            material()->setSpecularColor(QColor::fromRgbF(0.3, 0.3, 0.3));
        }
        else{
            material()->setSpecularColor(Qt::transparent);
        }

        if (d->colorMode != SolidColor) {
            setCalculated(false);
        }
    }
}

void GraphicsPymolSurfaceItem::setCalculated(bool calculated)
{
    if (!calculated) {
        delete d->buffer;
        d->buffer = 0;
        d->maxVdwCalculated = false;
    }
}

Real GraphicsPymolSurfaceItem::maxVdwRadius()
{
    if (!d->maxVdwCalculated) {
        if (!d->radii.empty())
            d->maxVdwRadius = *std::max_element(d->radii.begin(), d->radii.end());
        else
            d->maxVdwRadius = 0;
        d->maxVdwCalculated = true;
    }

    return d->maxVdwRadius;
}

} // end chemkit namespace
