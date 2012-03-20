/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "xtcfileformat.h"

#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>

#include <rpc/types.h>
#include <rpc/xdr.h>
#include <rpc/types.h>
#include <chemkit/vector3.h>
#include <chemkit/unitcell.h>
#include <chemkit/trajectory.h>
#include <chemkit/trajectoryfile.h>
#include <chemkit/trajectoryframe.h>
#include <chemkit/cartesiancoordinates.h>

#include "../../3rdparty/xdrf/xdrf.h"

XtcFileFormat::XtcFileFormat()
    : chemkit::TrajectoryFileFormat("xtc")
{
}

bool XtcFileFormat::read(std::istream &input, chemkit::TrajectoryFile *file)
{
    // read data into temporary file
    std::string tempFileName = (boost::filesystem::temp_directory_path() /
                                boost::filesystem::unique_path()).string();

    std::ofstream ostream(tempFileName.c_str());
    ostream << input.rdbuf();
    size_t dataSize = ostream.tellp();
    ostream.close();

    XDR xdrs;
    xdropen(&xdrs, tempFileName.c_str(), "r");

    boost::shared_ptr<chemkit::Trajectory> trajectory = boost::make_shared<chemkit::Trajectory>();

    while(xdr_getpos(&xdrs) < dataSize){
        // read magic (should be '1995')
        int magic = 0;
        xdr_int(&xdrs, &magic);
        if(magic != 1995){
            break;
        }

        // read atom count
        int atomCount = 0;
        xdr_int(&xdrs, &atomCount);

        // set trajectory size
        if(trajectory->isEmpty() || trajectory->size() < size_t(atomCount)){
            trajectory->resize(atomCount);
        }

        // create new frame
        chemkit::TrajectoryFrame *frame = trajectory->addFrame();

        // read frame number
        int frameNumber = 0;
        xdr_int(&xdrs, &frameNumber);

        // read time
        float time = 0;
        xdr_float(&xdrs, &time);

        // read unit cell
        float box[3][3];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                xdr_float(&xdrs, &box[i][j]);
            }
        }

        chemkit::Vector3 x(box[0][0], box[0][1], box[0][2]);
        chemkit::Vector3 y(box[1][0], box[1][1], box[1][2]);
        chemkit::Vector3 z(box[2][0], box[2][1], box[2][2]);

        frame->setUnitCell(new chemkit::UnitCell(x * 10, y * 10, z * 10));

        // read coordinates
        std::vector<float> coordinateData(3 * atomCount);
        float precision = 1000.0f;
        xdr3dfcoord(&xdrs, &coordinateData[0], &atomCount, &precision);

        for(int i = 0; i < atomCount; i++){
            // multiply each coordinate by 10 to convert
            // from nanometers to angstroms
            chemkit::Point3 position(coordinateData[i*3+0] * 10,
                                     coordinateData[i*3+1] * 10,
                                     coordinateData[i*3+2] * 10);

            frame->setPosition(i, position);
        }
    }

    xdrclose(&xdrs);

    // remove temp file
    boost::filesystem::remove(tempFileName);

    if(trajectory->isEmpty()){
        return false;
    }

    file->setTrajectory(trajectory);

    return true;
}
