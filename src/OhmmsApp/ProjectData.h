//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_PROJECTDATA_H__
#define OHHMS_PROJECTDATA_H__

#include "OhmmsData/OhmmsElementBase.h"
#include <vector>
#include <string>
#include <iostream>
using namespace std;

namespace OHMMS {

  /**class ProjectData
   *\brief Encapsulate data for a project
   *
   * Default: myName = "Project"
   * Should not modify the name, since composite types, such as MDRunData, use the name.
   *
   */
  struct ProjectData: public OhmmsElementBase {

    /// constructor
    ProjectData(const char* aname="project");
    
    bool get(ostream& os) const;
    bool put(istream& is);
    bool put(xmlNodePtr cur);
    void reset();
    
    ///increment a series number and reset m_projectroot
    void advance();
    
    ///returns the name of the project
    inline const char* CurrentRoot() const { return m_projectroot.c_str();}

    ///returns the name of the project
    inline const char* NextRoot() const { return m_nextroot.c_str();}

    ///\param oldroot is composed by the m_title and m_series
    bool PreviousRoot(string& oldroot) const;

    ///title of the project
    string m_title;

    ///user name
    string m_user;

    ///name of the host where the job is running
    string m_host;

    ///date when the job is executed
    string m_date;

    ///root for all the output engines
    string m_projectroot;

    ///root for the next run
    string m_nextroot;

    ///series index
    int m_series;

    ///the xml node for <Project/>
    xmlNodePtr m_cur;
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
