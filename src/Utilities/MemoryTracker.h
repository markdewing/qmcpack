//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef MEMORY_TRACKER_H
#define MEMORY_TRACKER_H

#include <malloc.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <set>


// Forward declarations
class Libxml2Document;
struct _xmlNode;
typedef _xmlNode* xmlNodePtr;

struct tracked_mem_object
{
  size_t num_objects;
  size_t num_objects_created;
  size_t current_bytes;
  size_t total_bytes;
  int first_create_time;
  int last_delete_time;
  tracked_mem_object (size_t size) :
    num_objects(1), num_objects_created(1),current_bytes(size),total_bytes(size),
    first_create_time(0),last_delete_time(0)
  {
  }
  tracked_mem_object() : num_objects(0), num_objects_created(0), current_bytes(0),
                         total_bytes(0), first_create_time(0), last_delete_time(0)
  { }
};


extern std::vector<std::string> thread_tags;
#pragma omp threadprivate(thread_tags)
  

class MemoryTrackerClass
{
private:
  std::map<std::string,tracked_mem_object> mem_map;
  std::map<void*,std::pair<std::string,size_t> > pointer_map;

  std::map<std::string,tracked_mem_object> tagged_mem_map;
  //std::set<std::string> tags;

  void get_tagged_name(const std::string &name, std::string &tagged_name);

public:
  void add(void *p, size_t bytes, const std::string &name);
  void resize(const std::string &name, long int delta_bytes, size_t new_bytes);

  void startTag(const std::string &tag);
  void endTag(const std::string &tag);


  void remove(void *p);
  void report();
  void output_memory(Libxml2Document &doc, xmlNodePtr root);
  
  void make_unique_name(std::string &name);
  void compress_memory_objects(std::map<std::string, tracked_mem_object> &compressed_obj_map);
  void fill_in_hierarchy(const std::map<std::string, tracked_mem_object> &compressed_obj_map,
      std::map<std::string, tracked_mem_object> &filled_map);


};

extern MemoryTrackerClass MemoryTracker;


#endif
