#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <bits/stdc++.h>

using namespace std;
using namespace Rcpp;

// simple implementation of a BK tree (https://en.wikipedia.org/wiki/BK-tree)
// for std::string elements and Hamming distance
// based on:
//   https://github.com/talhasaruhan/fuzzy-search/blob/master/bktree.hpp
class BKtree {
public:

  // tree node
  class node {
  public:
    int index; // items[index] in tree contains string of this node
    unordered_map<int, node*> children;

    // node constructor / destructor
    node(int i): index(i) {}
    ~node() {
      for (const pair<int, node*>&p : children) {
        delete p.second;
      }
    }
    
    // node capacity (memory consumption of node and its children)
    int capacity(){
      int sz = sizeof(index);
      sz += ((sizeof(int) + sizeof(node*)) * children.size()) + sizeof(children);
      for (pair<const int, node*>&child : children) {
        sz += child.second->capacity();
      }
      return sz;
    }
  };
  
  size_t size; // number of strings in tree
  
  // tree constructors
  BKtree(string str): size(1) { // from a string
    items.push_back(str);
    root = new node(0);
  }

  BKtree(): size(0){ // empty tree
    root = nullptr;
  }
  
  BKtree(const vector<string>& v) { // from a vector of strings
    items = v;
    size = v.size();
    _build_from_items(v);
  }

  // add a new string into to tree
  void insert(string str) {
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
  void remove(string str) {
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

  // search elements within distance tol
  vector<string> search(string str, int tol) {
    if (size <= 0) {
      // tree is empty, nothing to return
      return vector<string>(0);
    }
    vector<string> vec; // results vector
    _search(str, tol, root, vec); // start recursive search at the root
    return vec;
  }
  
  // does the tree contain a string within distance tol of str?
  bool has(string str, int tol = 0) {
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
  void print(node* r = nullptr, node* n = nullptr, int depth = 0) {
    Rcout << size << "\n";
    if (size > 0) {
      _print(root, root, 0); // recursive node printing starting from root at zero identation
    }
  }

  // tree capacity (memory consumption of tree plus root node and its children)    
  int capacity() {
    int sz = sizeof(size) + sizeof(node*) + sizeof(items) + sizeof(deleted);
    for (string &s: items) {
      sz += s.size() * sizeof(char);
    }
    for (const string &s: deleted) {
      sz += s.size() * sizeof(char);
    }
    if (size > 0) {
      sz += root->capacity();
    }
    return sz;
  }
  
  // calculate levenshtein distance between pair of strings
  static int levenshtein_distance(string str1, string str2){
    int m = str1.length(), n = str2.length();
    
    int** dn = new int*[m + 1];
    for(int i = 0; i < m + 1; ++i) {
      dn[i] = new int[n + 1];
    }
    
    for (int i = 0; i < m + 1; ++i) {
      for (int j = 0; j < n + 1; ++j) {
        if (i == 0) {
          dn[0][j] = j; // deletions at the start of str2

        } else if (j==0) {
          dn[i][0] = i; // deletions at the start of str1
        
        } else if (str1[i-1] == str2[j-1]) { // match of str1[i-1] and str2[j-1]
          dn[i][j] = dn[i-1][j-1];

        } else { // mismatch between str1[i-1] and str2[j-1] -> find minimal source
          dn[i][j] = 1 + min(
            dn[i-1][j], // deletion in str1
                   min(dn[i][j-1], // deletion in str2
                       dn[i-1][j-1]) // mismatch
          );
        }
      }
    }
    
    int d = dn[m][n];
    for (int i = 0; i < m + 1; ++i) {
      delete [] dn[i];
    }
    delete [] dn;
    
    return d;
  }
  
  // calculate hamming distance between pair of strings of equal lengths
  static int hamming_distance(string str1, string str2){
    int d = 0;
    
    for (size_t i = 0; i <= str1.length(); i++) {
      if (str1[i] != str2[i]) {
        d++;
      }
    }
    
    return d;
  }
  
  inline vector<string>::iterator begin() noexcept { return items.begin(); }
  inline vector<string>::const_iterator cbegin() const noexcept { return items.cbegin(); }
  inline vector<string>::iterator end() noexcept { return items.end(); }
  inline vector<string>::const_iterator cend() const noexcept { return items.end(); }
  
private:
  vector<string> items;
  node* root;
  unordered_set<string> deleted;
  
  // build new tree from string elements in vector v
  void _build_from_items(const vector<string>& v) {
    if (size < 1) {
      // empty vector, just create an empty tree
      root = nullptr;
      return;

    } else if (size == 1) {
      // single element -> root at index zero
      root = new node(0);
      return;

    } else {
      // multiple elements -> create rood and add rest sequentially
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
    items = newitems;
    
    // ... and empty deleted
    deleted.clear();
    
    // build up new tree from cleaned items
    delete root;
    size = items.size();
    _build_from_items(items);
  }
  
  // add items[index] as a new node to tree
  bool _add_from_items(int index) {
    string str = items[index];
    bool res = false;
    
    if(root == nullptr) { // tree is empty -> make it root
      root = new node(index);
      res = true;
      
    } else {
      node* t = root; // start with root
      node* new_node = new node(index); // create new node
      
      // distance of new item to current node
      // int dist = levenshtein_distance(items[t->index], str);
      int dist = hamming_distance(items[t->index], str);
      
      // while current node t already has a child at that distance dist
      while (t->children.find(dist) != t->children.end()) {
        // descend to the children of that child
        t = t->children.find(dist)->second;
        // ... and calculate new distance
        // dist = levenshtein_distance(items[t->index], str);
        dist = hamming_distance(items[t->index], str);
      }
      
      // t now points to a node without children at distance dist
      if (dist > 0) { // non-zero distance -> insert new string as child at distance dist
        t->children.insert(pair<int, node*>(dist, new_node));
        res = true;

      } else { // zero distance -> found new string already in the tree: don't insert and return false
        res = false;
      }
    }
    
    return res;
  }
  
  // recursively search for elements within tol of str in subtree of node r, add results to vec
  void _search(string str, int tol, node* r, vector<string>& vec) {
    string r_val = items[r->index];
    // int dist = levenshtein_distance(r_val, str);
    int dist = hamming_distance(r_val, str);
    
    if (dist <= tol) { // found element within tolerance
      // only add to results if not deleted
      if (deleted.find(r_val) == deleted.end()) {
        vec.push_back(r_val);
      }
    }
    
    // apply BK tree triangle inequality to
    // calculate boundaries for potential further hits relative to current hit
    int gte = dist - tol, lte = dist + tol;
    
    for (const pair<int, node*>& p: r->children) {
      if(p.first >= gte && p.first <= lte) {
        // ... and recurse search on children in that range
        _search(str, tol, p.second, vec);
      }
    }
  }
  
  // recursively check if subtree of node r has an element within tol of str
  bool _has(string str, int tol, node* r) {
    string r_val = items[r->index];
    // int dist = levenshtein_distance(r_val, str);
    int dist = hamming_distance(r_val, str);
    
    if (dist <= tol) { // found an element within tolerance
      // only return true if not deleted
      if (deleted.find(r_val) == deleted.end()) {
        return true;
      }
    }
    
    // apply BK tree triangle inequality to
    // calculate boundaries for potential further hits relative to current hit
    int gte = dist - tol, lte = dist + tol;
    
    for (const pair<int, node*>&p: r->children) {
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
    
    string n_val = items[n->index], r_val = items[r->index];
    if(deleted.find(n_val) != deleted.end()) {
      // if node n is deleted, prefix it with D*
      Rcout << "D*";
    }
    // Rcout << n_val << ", " << levenshtein_distance(r_val, n_val) << endl;
    Rcout << n_val << ", " << hamming_distance(r_val, n_val) << endl;
    for (const pair<int, node*>&x: n->children) {
      // recursively print all children of n
      _print(n, x.second, depth + 1);
    }
  }
};


// global instance of BK tree
BKtree tree;

// accessor functions exported to R
// get nuber of elements stored in tree
// [[Rcpp::export]]
size_t bk_size() {
  return tree.size;
}

// remove all elements from the tree
// [[Rcpp::export]]
size_t bk_clear() {
  tree.remove_all();
  return tree.size;
}

// create a new tree from a vector of elements
// [[Rcpp::export]]
size_t bk_new(std::vector<std::string> seqs, bool verbose = false) {
  tree.remove_all();
  tree = BKtree(seqs);
  if (verbose)
    Rcout << "new tree of size: " << tree.size << " (capacity: " << tree.capacity() << ")" << std::endl;
  return tree.size;
}

// add one or several elements to the tree
// [[Rcpp::export]]
size_t bk_add(std::vector<std::string> seqs, bool verbose = false) {
  size_t start = tree.size;
  for (size_t i = 0; i < seqs.size(); i++) {
    tree.insert(seqs[i]);
  }
  if (verbose)
    Rcout << "added " << (tree.size - start) << " elements, new size: " << tree.size << " (capacity: " << tree.capacity() << ")" << std::endl;
  return tree.size;
}

// remove one or several elements from the tree (trigger a re-build if needed)
// [[Rcpp::export]]
size_t bk_remove(std::vector<std::string> seqs, bool verbose = false) {
  size_t start = tree.size;
  for (size_t i = 0; i < seqs.size(); i++) {
    tree.remove(seqs[i]);
  }
  if (verbose)
    Rcout << "removed " << (start - tree.size) << " elements, new size: " << tree.size << " (capacity: " << tree.capacity() << ")" << std::endl;
  return tree.size;
}

// print tree to console (including elements deleted since last rebuild, flagged by "D*" prefix)
// [[Rcpp::export]]
void bk_print() {
  Rcout << "current size: " << tree.size << " (capacity: " << tree.capacity() << ")" << std::endl;
  tree.print();
}

// check if tree contains element seq
// [[Rcpp::export]]
bool bk_has(std::string seq, int tol = 0) {
  return tree.has(seq, tol);
}

// search elements within tol of seq
// [[Rcpp::export]]
std::vector<std::string> bk_search(std::string seq, int tol = 1, bool verbose = false) {
  std::vector<std::string> res = tree.search(seq, tol);
  if (verbose)
    Rcout << "found " << res.size() << " elements within distance " << tol << std::endl;
  return res;
}

// group elements by "greedy_clustering"
// [[Rcpp::export]]
List greedy_clustering(int tol = 1, bool verbose = false) {
  List L;
  std::string qseq;
  std::vector<std::string> gseqs;
  
  // start querying in the order of tree.items (assume ordered in a meaningful way)
  double start = (double)tree.size;
  while (tree.size > 0) {
    qseq = tree.first();
    gseqs = tree.search(qseq, tol);
    L.push_back(gseqs, qseq);
    bk_remove(gseqs, false);
  }

  return L;
}