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

#ifndef ROTATABLEBONDSPLUGIN_H
#define ROTATABLEBONDSPLUGIN_H

#include <chemkit/plugin.h>
#include <chemkit/moleculardescriptor.h>

class RotatableBondsPlugin : public chemkit::Plugin
{
    Q_OBJECT

    public:
        RotatableBondsPlugin();
        ~RotatableBondsPlugin();

        static chemkit::MolecularDescriptor* createRotableBondsDescriptor();
};

#endif // ROTATABLEBONDSPLUGIN_H