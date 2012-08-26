/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include <chemkit/plugin.h>
#include <chemkit/moleculardescriptor.h>
#include <chemkit/forcefieldenergydescriptor.h>

#include "amberatomtyper.h"
#include "amberforcefield.h"

#ifdef CHEMKIT_WITH_MD_IO
#include "mdcrdfileformat.h"
#endif

class AmberPlugin : public chemkit::Plugin
{
public:
    AmberPlugin()
        : chemkit::Plugin("amber")
    {
        CHEMKIT_REGISTER_ATOM_TYPER("amber", AmberAtomTyper);
        CHEMKIT_REGISTER_FORCE_FIELD("amber", AmberForceField);
        registerPluginClass<chemkit::MolecularDescriptor>("amber-energy", createAmberEnergyDescriptor);

        #ifdef CHEMKIT_WITH_MD_IO
        CHEMKIT_REGISTER_TRAJECTORY_FILE_FORMAT("mdcrd", MdcrdFileFormat);
        CHEMKIT_REGISTER_TRAJECTORY_FILE_FORMAT("trj", MdcrdFileFormat);
        #endif
    }

    static chemkit::MolecularDescriptor* createAmberEnergyDescriptor()
    {
        return new chemkit::ForceFieldEnergyDescriptor<AmberForceField>("amber-energy");
    }
};

CHEMKIT_EXPORT_PLUGIN(amber, AmberPlugin)
