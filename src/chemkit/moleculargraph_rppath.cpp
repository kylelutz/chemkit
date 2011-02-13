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

// This file implements the RP-Path ring perception algorithm. See
// [Lee 2009].

#include "molecule.h"
#include "moleculargraph.h"

namespace chemkit {

namespace {

class DistanceMatrix
{
    public:
        DistanceMatrix(int size);
        ~DistanceMatrix();

        int operator()(int i, int j) const;
        int& operator()(int i, int j);

    private:
        int m_size;
        int *m_values;
};

DistanceMatrix::DistanceMatrix(int size)
{
    m_size = size;
    m_values = new int[size*size];
    qMemSet(m_values, 0, size*size*sizeof(*m_values));
}

DistanceMatrix::~DistanceMatrix()
{
    delete [] m_values;
}

int DistanceMatrix::operator()(int i, int j) const
{
    return m_values[(i*m_size)+j];
}

int& DistanceMatrix::operator()(int i, int j)
{
    return m_values[(i*m_size)+j];
}

// path-included distance matrix
class PidMatrix
{
    public:
        PidMatrix(int size);
        PidMatrix(const Molecule *molecule);
        ~PidMatrix();

        void setPath(int i, int j, QVector<int> path);
        void setPath(int i, int j, QVector<QVector<int> > paths);
        void appendPath(int i, int j, QVector<QVector<int> > paths);
        void clearPath(int i, int j);

        const Atom* atom(int index) const;
        void setDistance(int i, int j, int distance);
        int distance(int i, int j) const;

        void setShortPath(int i, int j, const QVector<int> &path);
        void setShortPaths(int i, int j, QVector<QVector<int> > paths);
        void appendShortPath(int i, int j, const QVector<int> &path);
        void setLongPath(int i, int j, const QVector<int> &path);
        void appendLongPath(int i, int j, const QVector<int> &path);

        QVector<QVector<int> > paths(int i, int j) const { return m_values[(i*m_size)+j]; }

        QVector<QVector<int> > operator()(int i, int j) const;
        QVector<QVector<int> >& operator()(int i, int j);

        QVector<QVector<int> > splice(int i, int j, int k);

    private:
        int m_size;
        QVector<QVector<int> > *m_values;
};

PidMatrix::PidMatrix(int size)
{
    m_size = size;
    m_values = new QVector<QVector<int> >[size*size];
}

PidMatrix::~PidMatrix()
{
    delete [] m_values;
}

void PidMatrix::setPath(int i, int j, QVector<QVector<int> > paths)
{
    m_values[(i*m_size)+j].clear();
    m_values[(i*m_size)+j] = paths;
}

void PidMatrix::setPath(int i, int j, QVector<int> path)
{
    m_values[(i*m_size)+j].clear();
    m_values[(i*m_size)+j].append(path);
}

void PidMatrix::appendPath(int i, int j, QVector<QVector<int> > paths)
{
    foreach(QVector<int> path, paths){
        m_values[(i*m_size)+j].append(path);
    }
}

QVector<QVector<int> > PidMatrix::operator()(int i, int j) const
{
    return m_values[(i*m_size)+j];
}

QVector<QVector<int> >& PidMatrix::operator()(int i, int j)
{
    return m_values[(i*m_size)+j];
}

QVector<QVector<int> > PidMatrix::splice(int i, int j, int k)
{
    QVector<QVector<int> > splicedPaths;

    QVector<QVector<int> > ij_paths = paths(i, j);
    QVector<QVector<int> > jk_paths = paths(j, k);

    if(ij_paths.isEmpty() && jk_paths.isEmpty()){
        QVector<int> path;
        path.append(j);
        splicedPaths.append(path);
    }
    else if(ij_paths.isEmpty()){
        foreach(QVector<int> jk_path, jk_paths){
            QVector<int> path;
            path.append(j);
            path += jk_path;
            splicedPaths.append(path);
        }
    }
    else if(jk_paths.isEmpty()){
        foreach(QVector<int> ij_path, ij_paths){
            QVector<int> path = ij_path;
            path.append(j);
            splicedPaths.append(path);
        }
    }
    else{
        foreach(QVector<int> ij_path, ij_paths){
            foreach(QVector<int> jk_path, jk_paths){
                QVector<int> path = ij_path + (QVector<int>()<<j) + jk_path;
                splicedPaths.append(path);
            }
        }
    }

    return splicedPaths;
}

class RingCandidate
{
    public:
        RingCandidate(int size, int start, int end) : m_size(size), m_start(start), m_end(end) { }

        int size() const { return m_size; }
        int start() const { return m_start; }
        int end() const { return m_end; }

        static bool compareSize(RingCandidate *a, RingCandidate *b) { return a->size() < b->size(); }

    private:
        int m_size;
        int m_start;
        int m_end;
};

bool isUnique(const QVector<int> &path, const QList<QVector<int> > &sssr)
{

    QSet<int> pathSet = QSet<int>::fromList(path.toList());

    QVector<int> ringPath;
    foreach(ringPath, sssr){
        QSet<int> ringPathSet = QSet<int>::fromList(ringPath.toList());
        if(ringPathSet.subtract(pathSet).isEmpty()){
            return false;
        }
    }

    QSet<QPair<int, int> > ringBonds;
    foreach(ringPath, sssr){
        for(int i = 0; i < ringPath.size()-1; i++){
            ringBonds.insert(qMakePair(qMin(ringPath[i], ringPath[i+1]),
                                       qMax(ringPath[i], ringPath[i+1])));
        }
        ringBonds.insert(qMakePair(qMin(ringPath.first(), ringPath.last()),
                                   qMax(ringPath.first(), ringPath.last())));
    }

    int uniqueBonds = 0;
    for(int i = 0; i < path.size()-1; i++){
        if(!ringBonds.contains(qMakePair(qMin(path[i], path[i+1]),
                                         qMax(path[i], path[i+1])))){
            uniqueBonds++;
        }
    }
    if(!ringBonds.contains(qMakePair(qMin(path.first(), path.last()),
                                     qMax(path.first(), path.last())))){
        uniqueBonds++;
    }

    if(uniqueBonds == 0){
        return false;
    }

    return true;
}

} // end anonymous namespace

// The sssr_rpPath() method returns the smallest set of smallest rings in a
// molecular graph using the RP-Path algorithm. The graph is expected to
// contain a single fragment and to have all terminal nodes removed (i.e. all
// verticies should have degree >= 2).
QList<Ring *> MolecularGraph::sssr_rpPath(const MolecularGraph *graph)
{
    int n = graph->size();

    int ringCount = graph->bondCount() - graph->atomCount() + 1;
    if(ringCount < 1)
        return QList<Ring *>();

    // algorithm 1 - create the distance and pid matricies
    DistanceMatrix D(n);
    PidMatrix P(n);
    PidMatrix Pt(n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i == j){
                D(i, j) = 0;
            }
            else if(graph->adjacent(i, j)){
                D(i, j) = 1;
            }
            else{
                D(i, j) = INT_MAX/2; // ~ infinity
            }
        }
    }

    for(int k = 0; k < n; k++){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == j || i == k || k == j)
                    continue;

                if(D(i, j) > D(i, k) + D(k, j)){
                    if(D(i, j) == D(i, k) + D(k, j) + 1){
                        Pt(i, j) = P(i, j);
                    }
                    else{
                        Pt(i, j).clear();
                    }

                    D(i, j) = D(i, k) + D(k, j);
                    P(i, j) = P.splice(i, k, j);
                }
                else if(D(i, j) == D(i, k) + D(k, j)){
                    P(i, j) += P.splice(i, k, j);
                }
                else if(D(i, j) == D(i, k) + D(k, j) - 1){
                    Pt(i, j) += P.splice(i, k, j);
                }
            }
        }
    }

    // algorithm 2 - create the ring candidate set
    QList<RingCandidate *> ringCandidates;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            if(P(i, j).size() == 1 && Pt(i, j).size() == 0){
                continue;
            }
            else{
                int size;

                if(P(i, j).size() > 1){
                    size = 2 * D(i, j);
                }
                else{
                    size = 2 * D(i, j) + 1;
                }

                if(size > 2){
                    ringCandidates.append(new RingCandidate(size, i, j));
                }
            }
        }
    }
    qSort(ringCandidates.begin(), ringCandidates.end(), &RingCandidate::compareSize);

    // algorithm 3 - find sssr from the ring candidate set
    QList<QVector<int> > sssr;
    foreach(RingCandidate *ringCandidate, ringCandidates){

        // odd sized ring
        if(ringCandidate->size() & 1){
            for(int j = 0; j < Pt(ringCandidate->start(), ringCandidate->end()).size(); j++){
                QVector<int> path;
                path.append(ringCandidate->start());
                path += Pt(ringCandidate->start(), ringCandidate->end())[j];
                path.append(ringCandidate->end());
                if(P(ringCandidate->end(), ringCandidate->start()).size()){
                    path += P(ringCandidate->end(), ringCandidate->start())[0];
                }

                // check if ring is valid
                bool valid = true;
                foreach(int atom, path){
                    if(path.count(atom) > 1){
                        valid = false;
                        break;
                    }
                }

                if(valid && isUnique(path, sssr)){
                    sssr.append(path);
                    break;
                }
            }
        }

        // even sized ring
        else{
            for(int j = 0; j < P(ringCandidate->start(), ringCandidate->end()).size()-1; j++){
                QVector<int> path;
                path.append(ringCandidate->start());
                path += P(ringCandidate->start(), ringCandidate->end())[j];
                path.append(ringCandidate->end());
                path += P(ringCandidate->end(), ringCandidate->start())[j+1];

                // check if ring is valid
                bool valid = true;
                foreach(int atom, path){
                    if(path.count(atom) > 1){
                        valid = false;
                         break;
                    }
                }

                if(valid && isUnique(path, sssr)){
                    sssr.append(path);
                    break;
                }
            }
        }

        if(sssr.size() == ringCount){
            break;
        }
    }

    qDeleteAll(ringCandidates);

    QList<Ring *> rings;
    foreach(const QVector<int> &path, sssr){
        QList<Atom *> atomPath;

        foreach(int atomIndex, path){
            atomPath.append(const_cast<Atom *>(graph->atom(atomIndex)));
        }

        rings.append(new Ring(atomPath));
    }

    return rings;
}

} // end chemkit namespace
