//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <cassert>
#include <Utilities/MemoryTracker.h>
#include <Utilities/MemoryTracker_C.h>
#include <OhmmsData/Libxml2Doc.h>
#include <sstream>
#include "config.h"

int MemoryTrackerUniqueId=0;


MemoryTrackerClass MemoryTracker;

std::vector<std::string> thread_tags;

int alloc_clock = 1;

void MemoryTrackerClass::add(void *p, size_t bytes, const std::string &name)
{
#ifdef ENABLE_MEMORY_TRACKING
  if (active) 
  {
    #pragma omp critical
    {
      pointer_map[p] = std::pair<std::string,size_t> (name,bytes);
      std::map<std::string, tracked_mem_object>::iterator iter = mem_map.find(name);
      if (iter == mem_map.end())
      {
        tracked_mem_object obj(bytes);
        obj.first_create_time = alloc_clock;
        mem_map[name] = obj;
      }
      else
      {
        tracked_mem_object &obj = iter->second;
        obj.num_objects++;
        obj.num_objects_created++;
        obj.total_bytes += bytes;
        obj.current_bytes += bytes;
      }

      std::string tagged_name;
      get_tagged_name(name, tagged_name);
      //printf("adding %lu %s  tagged %s\n",bytes,name.c_str(), tagged_name.c_str());
      std::map<std::string, tracked_mem_object>::iterator tagged_iter = tagged_mem_map.find(tagged_name);
      if (tagged_iter == tagged_mem_map.end())
      {
        tracked_mem_object obj(bytes);
        obj.first_create_time = alloc_clock;
        alloc_clock++;
        tagged_mem_map[tagged_name] = obj;
      }
      else
      {
        tracked_mem_object &obj = tagged_iter->second;
        obj.num_objects++;
        obj.num_objects_created++;
        obj.total_bytes += bytes;
        obj.current_bytes += bytes;
      }
    }
  }
#endif
}

void MemoryTrackerClass::resize(const std::string &name, long int delta_bytes, size_t new_bytes)
{
#ifdef ENABLE_MEMORY_TRACKING
  if (active) 
  {
    #pragma omp critical
    {
      std::map<std::string, tracked_mem_object>::iterator iter = mem_map.find(name);
      if (iter == mem_map.end())
      {
        tracked_mem_object obj(delta_bytes);
        obj.first_create_time = alloc_clock;
        mem_map[name] = obj;
      }
      else
      {
        tracked_mem_object &obj = iter->second;
        obj.total_bytes += delta_bytes;
        obj.current_bytes += delta_bytes;
        if (obj.current_bytes != new_bytes) {
          printf("  error in resize: new bytes = %ld  computed bytes = %ld\n",new_bytes,obj.total_bytes);
        }
      }

      std::string tagged_name;
      get_tagged_name(name, tagged_name);
      //printf("resizing (by %ld) %s  tagged %s new bytes%ld\n",delta_bytes, name.c_str(), tagged_name.c_str(),new_bytes);
      std::map<std::string, tracked_mem_object>::iterator tagged_iter = tagged_mem_map.find(tagged_name);
      if (tagged_iter == tagged_mem_map.end())
      {
        tracked_mem_object obj(delta_bytes);
        obj.first_create_time = alloc_clock;
        alloc_clock++;
        tagged_mem_map[tagged_name] = obj;
      }
      else
      {
        tracked_mem_object &obj = tagged_iter->second;
        obj.total_bytes += delta_bytes;
        obj.current_bytes += delta_bytes;
      }
    }
  }
#endif
}

void MemoryTrackerClass::remove(void *p)
{
  remove(p, "", 0);
}

void MemoryTrackerClass::remove(const std::string &name, size_t bytes)
{
  remove(NULL, name, bytes);
}

void MemoryTrackerClass::remove(void *p, const std::string &name, size_t bytes)
{
#ifdef ENABLE_MEMORY_TRACKING
  if (active) 
  {
    #pragma omp critical
    {
      std::string alloc_name = name;
      size_t alloc_bytes = bytes;
      if (p != NULL)
      {
        std::map<void*,std::pair<std::string,size_t> >::iterator piter;
        piter = pointer_map.find(p);
        if (piter == pointer_map.end())
        {
          fprintf (stderr, "Attempt to free a pointer not in the map.\n");
          abort();
        }
        alloc_name = piter->second.first;
        alloc_bytes = piter->second.second;
        // fprintf (stderr, "Deallocating %s from GPU memory of size %ld at pointer %p.\n",
        //         name.c_str(), bytes, p);
        pointer_map.erase(piter);
      }

      std::map<std::string, tracked_mem_object>::iterator iter = mem_map.find(alloc_name);
  #if 0
      if (iter == mem_map.end())
      {
        fprintf (stderr, "Error:  Memory object %s not found "
                 "in the memory map\n", name.c_str());
        abort();
      }
  #endif
      if (iter != mem_map.end())
      {
        tracked_mem_object &obj = iter->second;
        obj.last_delete_time = alloc_clock;
        if (obj.num_objects == 1)
        {
          assert (alloc_bytes == obj.current_bytes);
        }
        obj.num_objects--;
        obj.current_bytes -= alloc_bytes;
      }

  #if 1
      std::string tagged_name;
      get_tagged_name(alloc_name, tagged_name);
      //printf("Removing %s  tagged %s\n",name.c_str(), tagged_name.c_str());
      std::map<std::string, tracked_mem_object>::iterator tagged_iter = tagged_mem_map.find(tagged_name);
      if (tagged_iter != tagged_mem_map.end())
      {
        tracked_mem_object &obj = tagged_iter->second;
        obj.last_delete_time = alloc_clock;
        alloc_clock++;
        if (obj.num_objects == 1)
        {
          assert (alloc_bytes == obj.current_bytes);
        }
        obj.num_objects--;
        obj.current_bytes -= alloc_bytes;
      }
  #endif
    }
  }
#endif
}

void
MemoryTrackerClass::startTag(const std::string &tag)
{
#ifdef ENABLE_MEMORY_TRACKING
  if (active)
  {
    //#pragma omp critical
    //thread_tags.insert(tag);
    thread_tags.push_back(tag);
  }
#endif
}

void
MemoryTrackerClass::endTag(const std::string &tag)
{
#ifdef ENABLE_MEMORY_TRACKING
  if (active)
  {
    //#pragma omp critical
    //thread_tags.erase(tag);
    if (thread_tags.back() != tag) {
      printf("Error in MemoryTracker: Popping tag %s but should be %s\n",
             thread_tags.back().c_str(), tag.c_str());
    }
    thread_tags.pop_back();
  }
#endif
}

void
MemoryTrackerClass::get_tagged_name(const std::string &name, std::string &tagged_name)
{
  tagged_name = name;
  //std::set<std::string>::iterator it = thread_tags.begin();
  std::vector<std::string>::iterator it = thread_tags.begin();
  for (; it != thread_tags.end(); ++it)
  {
    //tagged_name += "," + *it;
    tagged_name += "/" + *it;
  }
}

void
MemoryTrackerClass::make_unique_name(std::string &name)
{
  #pragma omp critical
  {
    std::ostringstream tmp;
    tmp << name << MemoryTrackerUniqueId;
    name = tmp.str();
    MemoryTrackerUniqueId++;
  }
}


int
get_level(const std::string &stack_name)
{
  int level = 0;
  for (int i = 0; i < stack_name.length(); i++)
  {
    if (stack_name[i] == '/')
    {
      level++;
    }
  }
  return level;
}


std::string
get_leaf_name(const std::string &stack_name)
{
  int pos = stack_name.find_last_of('/');
  if (pos == std::string::npos)
  {
    return stack_name;
  }

  return stack_name.substr(pos+1, stack_name.length()-pos);
}


void pad_string(const std::string &in, std::string &out, int field_len)
{
    int len = in.size();
    int pad_len = std::max(field_len - len, 0);
    std::string pad_str(pad_len, ' ');
    out = in + pad_str;
}

void
MemoryTrackerClass::report()
{
#ifdef ENABLE_MEMORY_TRACKING
  fprintf (stdout, "Object name                                              Num objects   Current objects     Total bytes current bytes\n");
  std::map<std::string, tracked_mem_object>::iterator iter;
  size_t total_bytes = 0, total_num=0;
  for (iter=mem_map.begin(); iter != mem_map.end(); iter++)
  {
    fprintf (stdout, "%40s %8ld       %13ld %13ld %13ld   %4d %4d\n",
             iter->first.c_str(),
             iter->second.num_objects_created,
             iter->second.num_objects,
             iter->second.total_bytes,
             iter->second.current_bytes,
             iter->second.first_create_time,
             iter->second.last_delete_time);
    total_bytes += iter->second.total_bytes;
    total_num += iter->second.num_objects;
  }
#if 0
  fprintf (stdout, "%60s %8ld            %13ld\n", "Total",
           total_num, total_bytes);

  std::map<std::string, tracked_mem_object>::iterator tagged_iter;
  fprintf (stdout, "\nTagged Object name                                       Num objects            Total bytes\n");
  for (tagged_iter=tagged_mem_map.begin(); tagged_iter != tagged_mem_map.end(); tagged_iter++)
  {
    fprintf (stdout, "%40s %8ld      %13ld %13ld %13ld %4d %4d\n",
             tagged_iter->first.c_str(),
             tagged_iter->second.num_objects_created,
             tagged_iter->second.num_objects,
             tagged_iter->second.total_bytes,
             tagged_iter->second.current_bytes,
             tagged_iter->second.first_create_time,
             tagged_iter->second.last_delete_time);
  }
#endif



  std::map<std::string, tracked_mem_object> compressed_map;
  compress_memory_objects(compressed_map);
  std::map<std::string, tracked_mem_object> filled_map;
  fill_in_hierarchy(compressed_map, filled_map);
  std::map<std::string, tracked_mem_object>::iterator c_iter;


  fprintf (stdout, "\nTagged Object name                                       Num objects            Total bytes\n");
  for (c_iter=filled_map.begin(); c_iter != filled_map.end(); c_iter++)
  {
    const std::string &stack_name  = c_iter->first.c_str();
    int level = get_level(stack_name);
    const std::string &name = get_leaf_name(stack_name);
    std::string indent_str(2*level, ' ');
    std::string indented_str = indent_str + name;
    std::string padded_name_str;
    pad_string(indented_str, padded_name_str, 40);
    fprintf (stdout, "%40s %8ld      %13ld %13ld %13ld %4d %4d\n",
             padded_name_str.c_str(),
             c_iter->second.num_objects_created,
             c_iter->second.num_objects,
             c_iter->second.total_bytes,
             c_iter->second.current_bytes,
             c_iter->second.first_create_time,
             c_iter->second.last_delete_time);
  }
#endif
}

// Matches .*_[0-9]+, and returns the base value if matched
bool matches_object_suffix(const std::string &in, std::string &base)
{
  base = in;
  size_t u_pos = in.find_last_of('_');
  if (u_pos == std::string::npos || u_pos == in.length()-1)
  {
    return false;
  }

  for (size_t i = u_pos+1; i < in.length(); i++)
  {
    if (!std::isdigit(in[i]))
    {
      return false;
    }
  }
  base = in.substr(0, u_pos);
  return true;
}

void split_stack(const std::string &in, std::string &top, std::string &rest)
{
  size_t pos = in.find('/');
  top = in;
  if (pos != std::string::npos)
  {
    top = in.substr(0, pos);
    rest = in.substr(pos+1,rest.length()-pos);
  }
}



void
MemoryTrackerClass::compress_memory_objects(std::map<std::string, tracked_mem_object> &compressed_obj_map)
{
  std::map<std::string, tracked_mem_object>::iterator tagged_iter;
  for (tagged_iter=tagged_mem_map.begin(); tagged_iter != tagged_mem_map.end(); tagged_iter++)
  {
     // get names in stack  - strip the _# suffix
    const std::string &stack = tagged_iter->first;
    std::string top_name;
    std::string rest;
    split_stack(stack, top_name, rest);
    std::string base_name;
    matches_object_suffix(top_name, base_name);
    std::string compressed_stack = base_name + "/" + rest;
     
    const std::map<std::string, tracked_mem_object>::iterator &c_iter = compressed_obj_map.find(compressed_stack);
    if (c_iter == compressed_obj_map.end())
    {
      compressed_obj_map[compressed_stack] = tagged_iter->second;
    }
    else
    {
      tracked_mem_object &obj = tagged_iter->second;
      tracked_mem_object &cobj = c_iter->second;
      cobj.num_objects_created += obj.num_objects_created;
      cobj.num_objects += obj.num_objects;
      cobj.current_bytes += obj.current_bytes;
      cobj.total_bytes += obj.total_bytes;
      cobj.first_create_time = std::min(cobj.first_create_time, obj.first_create_time);
      cobj.last_delete_time = std::max(cobj.last_delete_time, obj.last_delete_time);
      compressed_obj_map[compressed_stack] = cobj;
    }
  }
}

void
MemoryTrackerClass::fill_in_hierarchy(const std::map<std::string, tracked_mem_object> &compressed_map,
      std::map<std::string, tracked_mem_object> &filled_map)
{
  std::map<std::string, tracked_mem_object>::const_iterator iter;
  for (iter=compressed_map.begin(); iter != compressed_map.end(); iter++)
  {
    std::string name = iter->first;
    const tracked_mem_object &obj = iter->second;

    int level = get_level(name);
    std::string stack_name;
    for (int i = 0; i < level+1; i++)
    {
      std::string top_name;
      std::string rest;
      split_stack(name, top_name, rest);
      if (stack_name == "")
      {
        stack_name = top_name;
      }
      else
      {
        stack_name = stack_name + '/' + top_name;
      } 
      name = rest;

      const std::map<std::string, tracked_mem_object>::iterator &iter = filled_map.find(stack_name);
      if (iter == filled_map.end())
      {
        filled_map[stack_name] = obj;
      }
      else
      {
        tracked_mem_object &cobj = iter->second;
        cobj.num_objects_created = std::max(obj.num_objects_created, cobj.num_objects_created);
        cobj.num_objects = std::max(obj.num_objects, cobj.num_objects);
        cobj.current_bytes += obj.current_bytes;
        cobj.total_bytes += obj.total_bytes;
        cobj.first_create_time = std::min(cobj.first_create_time, obj.first_create_time);
        cobj.last_delete_time = std::max(cobj.last_delete_time, obj.last_delete_time);
        filled_map[stack_name] = cobj;
      }
    }
  }
}

void
MemoryTrackerClass::output_memory(Libxml2Document &doc, xmlNodePtr root)
{
#ifdef ENABLE_MEMORY_TRACKING
  std::map<std::string, tracked_mem_object> compressed_map;
  compress_memory_objects(compressed_map);

  std::map<std::string, tracked_mem_object> filled_map;
  fill_in_hierarchy(compressed_map, filled_map);

  std::vector<std::string> names;
  std::map<std::string, tracked_mem_object>::const_iterator iter;
  for (iter=filled_map.begin(); iter != filled_map.end(); iter++)
  {
    names.push_back(iter->first);
  }

  xmlNodePtr memory_root = doc.addChild(root, "memory");

  std::vector<xmlNodePtr> node_stack;
  node_stack.push_back(memory_root);
  xmlNodePtr current_root = memory_root;

#if 1
  for (int i = 0; i < names.size(); i++)
  {
    const std::string &stack_name  = names[i];
    int level = get_level(stack_name)+1;
    const std::string &name = get_leaf_name(stack_name);

    xmlNodePtr mem_obj = doc.addChild(current_root, "object");
    tracked_mem_object &obj = filled_map[stack_name];
    doc.addChild(mem_obj, "name", name);
    doc.addChild(mem_obj, "number", obj.num_objects_created);
    doc.addChild(mem_obj, "current_number", obj.num_objects);
    doc.addChild(mem_obj, "bytes", obj.total_bytes);
    doc.addChild(mem_obj, "current_bytes", obj.current_bytes);
    doc.addChild(mem_obj, "first_create_time", obj.first_create_time);
    doc.addChild(mem_obj, "last_delete_time", obj.last_delete_time);

    int next_level = level;
    if (i+1 < names.size())
    {
      next_level = get_level(names[i+1])+1;
    }
    if (next_level > level)
    {
      xmlNodePtr next_node = doc.addChild(mem_obj, "includes");
      node_stack.push_back(next_node);
      current_root = next_node;
    }
    if (next_level < level)
    {
      for (int j = 0; j < level-next_level; j++)
      {
        node_stack.pop_back();
        current_root = node_stack.back();
      }
    }
  }
#endif
#endif
}


void MemoryTracker_add(void *ptr, size_t bytes, const char *name)
{
  MemoryTracker.add(ptr, bytes, name);
}

void MemoryTracker_remove(void *ptr)
{
  MemoryTracker.remove(ptr);
}
