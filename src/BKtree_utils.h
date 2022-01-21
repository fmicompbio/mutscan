#ifndef BKTREE_UTILS_HPP
#define BKTREE_UTILS_HPP

#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include "stringdist.h"

using namespace Rcpp;


// simple implementation of a BK tree (https://en.wikipedia.org/wiki/BK-tree)
// for std::string elements and Hamming (+ Shifts) or Levenshtein distance
// based on:
//   https://github.com/talhasaruhan/fuzzy-search/blob/master/bktree.hpp
class BKtree {
public:

  // tree node
  class node {
  public:
    int index; // items[index] in tree contains string of this node
    std::unordered_map<int, node*> children;

    // node constructor / destructor
    node(int i): index(i) {}
    ~node() {
      for (const std::pair<int, node*>&p : children) {
        delete p.second;
      }
    }
    
    // node capacity (memory consumption of node and its children)
    int capacity(){
      int sz = sizeof(index);
      sz += ((sizeof(int) + sizeof(node*)) * children.size()) + sizeof(children);
      for (std::pair<const int, node*>&child : children) {
        sz += child.second->capacity();
      }
      return sz;
    }
  };
  
  size_t size; // number of strings in tree
  
  // tree constructors
  BKtree(const std::string dmetric = "hamming", int mshift = -1) { // empty tree
    size = 0;
    metric = dmetric;
    max_absolute_shift = mshift;
    _set_distance_function_pointer();
    root = nullptr;
  }
  
  BKtree(const std::vector<std::string>& v,
         const std::string dmetric = "hamming", int mshift = -1) { // from a vector of strings
    metric = dmetric;
    max_absolute_shift = mshift;
    _set_distance_function_pointer();
    items = v;
    size = v.size();
    _build_from_items(v);
  }

  // add a new string into to tree
  void insert(std::string str) {
    size++;
    
    // remove it from deleted map if necessary
    if (deleted.find(str) != deleted.end()){
      deleted.erase(str);
      return;
    }
    
    // add it to items and get its index (zero-based)
    int index = items.size();
    items.push_back(str);
    // ... and try to add it to the tree as a new child node
    if(!_add_from_items(index)) {
      // if new string str was already in the tree, remove it again from items and size
      size--;
      items.pop_back();
    }
  }
  
  // remove a string from the tree
  // - string is only remembered as deleted, but not yet removed from the tree
  // - if too many deleted strings exist, rebuild whole tree from scratch
  void remove(std::string str) {
    if(deleted.find(str) != deleted.end() || !has(str)) {
      // string is already deleted or not in tree, nothing to do
      return;
    }

    // add string to deleted map and decrease size
    deleted.insert(str);
    size--;
    if (deleted.size() >= size/2) {
      if (size > 0) {
        // deleted is now quite large -> rebuild tree and release memory
        _rebuild();
      } else {
        // tree is empty now
        delete root;
        root = nullptr;
        items.clear();
        deleted.clear();
      }
    }
  }
  
  // remove all elements from tree
  void remove_all() {
    if (size > 0) {
      delete root;
      root = nullptr;
      size = 0;
      items.clear();
      deleted.clear();
    }
  }
  
  // get all (non-deleted) elements from the tree
  std::vector<std::string> get_all() {
    std::vector<std::string> all;
    for (size_t i = 0; i < items.size(); i++) {
      if (deleted.find(items[i]) == deleted.end()) {
        // item[i] is not deleted -> copy to all
        all.push_back(items[i]);
      }
    }
    return all;
  }

  // search elements within distance tol
  std::vector<std::string> search(std::string str, int tol) {
    if (size <= 0) {
      // tree is empty, nothing to return
      return std::vector<std::string>(0);
    }
    std::vector<std::string> vec; // results vector
    _search(str, tol, root, vec); // start recursive search at the root
    return vec;
  }
  
  // does the tree contain a string within distance tol of str?
  bool has(std::string str, int tol = 0) {
    if (size <= 0) {
      // tree is empty, return false
      return false;
    }
    return _has(str, tol, root); // start  recursive search at the root
  }
  
  // get first non-deleted element (in the order of items)
  std::string first() {
    std::string res = "";
    if (size > 0) {
      for (size_t i = 0; i < size; i++) {
        if (deleted.find(items[i]) == deleted.end()) {
          res = items[i];
          break;
        }
      }
    }
    return res;
  }

  // print tree on stdout
  // void print(node* r = nullptr, node* n = nullptr, int depth = 0) {
  void print() {
    Rcout << size << "\n";
    if (size > 0) {
      _print(root, root, 0); // recursive node printing starting from root at zero identation
    }
  }

  // tree capacity (memory consumption of tree plus root node and its children)    
  int capacity() {
    int sz = sizeof(size) + sizeof(node*) + sizeof(items) + sizeof(deleted);
    for (std::string &s: items) {
      sz += s.size() * sizeof(char);
    }
    for (const std::string &s: deleted) {
      sz += s.size() * sizeof(char);
    }
    if (size > 0) {
      sz += root->capacity();
    }
    return sz;
  }
  
  // return the distance metric used by a tree instance
  std::string distance_metric() {
    return metric;
  }
  
  // return the maximal shift used in hamming_shift_distance()
  int maximal_absolute_shift() {
    return max_absolute_shift;
  }
  
  inline std::vector<std::string>::iterator begin() noexcept { return items.begin(); }
  inline std::vector<std::string>::const_iterator cbegin() const noexcept { return items.cbegin(); }
  inline std::vector<std::string>::iterator end() noexcept { return items.end(); }
  inline std::vector<std::string>::const_iterator cend() const noexcept { return items.end(); }
  
private:
  std::vector<std::string> items;            // all strings added to the tree in the order of addition
                                             // (including potentially deleted ones)
  node* root;                                // pointer to root node
  std::unordered_set<std::string> deleted;   // nodes deleted from the tree but not yet removed
  std::string metric;                        // name of the distance metric to use
  int (*distance)(const std::string&, const std::string&, int); // pointer to function to calcluate string distance
  int max_absolute_shift;                    // maximum shift (only used for metric="hamming_shift")

  // set the distance function pointer according to metric
  void _set_distance_function_pointer() {
    if (metric == "hamming") {
      distance = &hamming_distance;

    } else if (metric == "hamming_shift") {
      distance = &hamming_shift_distance;
      
    } else if (metric == "levenshtein") {
      distance = &levenshtein_distance;

    } else {
      stop("unknown distance metric '%s'", metric);
    }
  }
  
  // build new tree from string elements in vector v
  void _build_from_items(const std::vector<std::string>& v) {
    if (size < 1) {
      // empty vector, just create an empty tree
      root = nullptr;
      return;

    } else if (size == 1) {
      // single element -> root at index zero
      root = new node(0);
      return;

    } else {
      // multiple elements -> create root and add rest sequentially
      root = new node(0);
      for (size_t i = 1; i < size; ++i) {
        if (!_add_from_items(i)) { // attach item[i] as new child
          // item[i] was already in the tree -> put it to the end of items
          swap(items[i], items[size-1]);
          // and forget about it
          size--;
          i--;
        }
      }
      // truncate items to unique items
      items.resize(size);
    }
  }
  
  // rebuild tree from items and release elements stored in deleted
  void _rebuild() {
    // reorder items (move deleted items to the end)
    // remark: this implementation perturbs the original order in items
    /*
    size_t i = 0, j = items.size() - 1;
    while (i <= j) {
      if (deleted.find(items[i]) != deleted.end()) {
        // item[i] is deleted -> move it to the end of items
        swap(items[i], items[j--]);
      } else {
        i++;
      }
    }
    // truncate items to just the non-deleted elements
    items.resize(i);
    */
    
    // alternative implementation (preserving order in items)
    std::vector<std::string> newitems;
    for (size_t i = 0; i < items.size(); i++) {
      if (deleted.find(items[i]) == deleted.end()) {
        // item[i] is not deleted -> copy to newitems
        newitems.push_back(items[i]);
      }
    }
    items.clear();
    items.insert(items.begin(), newitems.begin(), newitems.end());
    // items = newitems;
    
    // ... and empty deleted
    deleted.clear();
    
    // build up new tree from cleaned items
    delete root;
    size = items.size();
    _build_from_items(items);
  }
  
  // add items[index] as a new node to tree
  bool _add_from_items(int index) {
    std::string str = items[index];
    bool res = false;
    
    if(root == nullptr) { // tree is empty -> make it root
      root = new node(index);
      res = true;
      
    } else {
      node* t = root; // start with root
      node* new_node = new node(index); // create new node
      
      // distance of new item to current node
      int dist = distance(items[t->index], str, max_absolute_shift);
      
      // while current node t already has a child at that distance dist
      while (t->children.find(dist) != t->children.end()) {
        // descend to the children of that child
        t = t->children.find(dist)->second;
        // ... and calculate new distance
        dist = distance(items[t->index], str, max_absolute_shift);
      }
      
      // t now points to a node without children at distance dist
      if (dist > 0) { // non-zero distance -> insert new string as child at distance dist
        t->children.insert(std::pair<int, node*>(dist, new_node));
        res = true;

      } else { // zero distance -> found new string already in the tree: don't insert and return false
        res = false;
      }
    }
    
    return res;
  }
  
  // recursively search for elements within tol of str in subtree of node r, add results to vec
  void _search(std::string str, int tol, node* r, std::vector<std::string>& vec) {
    std::string r_val = items[r->index];
    int dist = distance(r_val, str, max_absolute_shift);
    
    if (dist <= tol) { // found element within tolerance
      // only add to results if not deleted
      if (deleted.find(r_val) == deleted.end()) {
        vec.push_back(r_val);
      }
    }
    
    // apply BK tree triangle inequality to
    // calculate boundaries for potential further hits relative to current hit
    int gte = dist - tol, lte = dist + tol;
    
    for (const std::pair<int, node*>& p: r->children) {
      if(p.first >= gte && p.first <= lte) {
        // ... and recurse search on children in that range
        _search(str, tol, p.second, vec);
      }
    }
  }
  
  // recursively check if subtree of node r has an element within tol of str
  bool _has(std::string str, int tol, node* r) {
    std::string r_val = items[r->index];
    int dist = distance(r_val, str, max_absolute_shift);
    
    if (dist <= tol) { // found an element within tolerance
      // only return true if not deleted
      if (deleted.find(r_val) == deleted.end()) {
        return true;
      }
    }
    
    // apply BK tree triangle inequality to
    // calculate boundaries for potential further hits relative to current hit
    int gte = dist - tol, lte = dist + tol;
    
    for (const std::pair<int, node*>&p: r->children) {
      if (p.first >= gte && p.first <= lte) {
        // ... and start recursive search on children in that range
        if (_has(str, tol, p.second)) {
          return true;
        }
      }
    }
    
    return false;
  }
  
  // recursively print subtree underneath node n:
  // first print node n (child of node r), indented by depth tabs
  // and the recurse on all children of n
  void _print(node* r = nullptr, node* n = nullptr, int depth = 0) {
    if (n == nullptr) {
      n = root;
    }
    if (r == nullptr) {
      r = root;
    }
    
    for (int i = 0; i < depth; ++i) {
      Rcout << "\t";
    }
    
    std::string n_val = items[n->index], r_val = items[r->index];
    if(deleted.find(n_val) != deleted.end()) {
      // if node n is deleted, prefix it with D*
      Rcout << "D*";
    }
    Rcout << distance(r_val, n_val, max_absolute_shift) << ": " << n_val << std::endl;
    for (const std::pair<int, node*>&x: n->children) {
      // recursively print all children of n
      _print(n, x.second, depth + 1);
    }
  }
};



// exposing the BKtree class to R using Rcpp-modules
// currently only used for unit testing
RCPP_MODULE(mod_BKtree) {
  class_<BKtree>( "BKtree" )

  .constructor()
  .constructor<std::string>()
  .constructor<std::string, int>()
  .constructor<std::vector<std::string>>()
  .constructor<std::vector<std::string>, std::string>()
  .constructor<std::vector<std::string>, std::string, int>()

  .field_readonly( "size", &BKtree::size, "number of elements stored in tree" )

  .method( "insert", &BKtree::insert, "insert a new element into tree" )
  .method( "remove", &BKtree::remove, "remove an element from tree" )
  .method( "remove_all", &BKtree::remove_all, "remove all elements from tree" )
  .method( "get_all", &BKtree::get_all, "get all elements stored in tree" )
  .method( "print", &BKtree::print, "print tree on console" )
  .method( "search", &BKtree::search, "search for elements within distance tolerance" )
  .method( "has", &BKtree::has, "check if there are elements within distance toelrance" )
  .method( "first", &BKtree::first, "get first (non-deleted) element from tree" )
  .method( "capacity", &BKtree::capacity, "calculate tree memory usage" )
  .method( "distance_metric", &BKtree::distance_metric, "return distance metric used by tree")
  .method( "maximal_absolute_shift", &BKtree::maximal_absolute_shift, "return maximal shift for metric='hamming_shift'")

  ;
}

// TODO: Attempting to export this function here gives a warning: 
// Invalid parameter: 'std::string&' for Rcpp::export attribute at BKtree_utils.hpp:467 
// For now, don't export it. 
int findClosestRefSeqFaster(std::string&, BKtree&, std::map<std::string, int>&, int&);

#endif
