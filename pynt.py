#!/usr/bin/python
import sys
import re

# Constants

TOKEN_NONE = 0
TOKEN_OPAR = 1
TOKEN_CPAR = 2
TOKEN_COLON = 3
TOKEN_COMMA = 4
TOKEN_SEMICOLON = 5
TOKEN_STRING = 6

# State transitions
#
#             NONE   OPAR   CPAR   COLON   COMMA   SEMICOLON   STRING
#  NONE       False  True   False  False   False   False       True
#  OPAR       False  True   False  False   False   False       True
#  CPAR       False  False  True   True    True    True        True
#  COLON      False  False  False  False   False   False       True
#  COMMA      False  True   False  False   False   False       True
#  SEMICOLON  False  False  False  False   False   False       False
#  STRING     False  False  True   True    True    True        False

lookup = [[False, True,  False, False, False, False, True ],
          [False, True,  False, False, False, False, True ],
          [False, False, True,  True,  True,  True,  True ],
          [False, False, False, False, False, False, True ],
          [False, True,  False, False, False, False, True ],
          [False, False, False, False, False, False, False],
          [False, False, True,  True,  True,  True,  False]]


def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

class CNode:
  '''Tree node structure'''

  def __init__(self,label = None):
    self.length = 0
    self.children = []
    self.parent = None
    self.label = None
    self.mark = False
    self.label = label

class CTree:
  '''Wrapper class for tree structure'''


  def __init__(self,newick):
    self.root = None
    self.leaves = []
    self.inner = []
    self.newick = None
    self.leaves_count = 0
    self.inner_count = 0

    # 'private' variables
    self._nwk = ""

    self.newick = stripwhite(newick)
    #print "Stripped newick: " + self.newick

    tokens = tokenize(newick)
    #print "Tokenized: " + str(tokens)
    rc = validate_trivial(tokens)
    if rc is True:
      rc,root = parse(tokens)

    if rc is False:
      #raise Exception('spam','eggs')
      raise ValueError('Incorrect formatted newick string',newick)
    else:
    #  print "Root children:"
      for i in root.children:
        #print "  " + str(i.label)
        continue
      self.root = root
      self._fill_internals(root)

  def _fill_internals(self,node):
    if len(node.children) == 0:
      self.leaves += [node]
      self.leaves_count += 1
      return

    for c in node.children:
      self._fill_internals(c)
    self.inner += [node]
    self.inner_count += 1

  def clear_mark(self):
    for node in self.leaves:
      node.mark=False
    for node in self.inner:
      node.mark=False

  def _xprint_tree(self,node):
    if len(node.children) == 0:
      self._nwk += node.label + ":" + str(node.length)
      return

    self._nwk += "("
    for c in node.children[:-1]:
      self._xprint_tree(c)
      self._nwk += ","
    self._xprint_tree(node.children[-1])
    self._nwk += ")"
    if node.label <> None:
      self._nwk += node.label
    self._nwk += ":" + str(node.length)

  def print_tree(self):
    self._nwk = ""
    self._xprint_tree(self.root)
    print self._nwk

  def get_node(self,label):
      
    for node in self.leaves:
      if node.label == label:
        return node
    for node in self.inner:
      if node.label == label:
        return node
    
    return None # In case the label in not in the tree

  def subtree(self, labels):
    
    self.clear_mark()
    
    # get leaf objects from labels
    leaf_list = []
    for label in labels:
      for node in self.leaves:
        if node.label == label:
          leaf_list += [node]
    
    for node in leaf_list:
      node.mark = True
      while node.parent != None:
        node = node.parent
        node.mark=True
     
    subroot = self.root
    
    while 1:
      count = 0
      for child in subroot.children:
        if child.mark == True:
          candidate = child
          count += 1
      if count != 1:
        break
      subroot = candidate
    
    self.clear_mark()

    return subroot


def parse(tokens):
  root = None
  parent = None
  node = None
  rc = True

  prev_token = TOKEN_NONE

  for i in range(len(tokens)):
    t = tokens[i]

    if t == '(':
      if not lookup[prev_token][TOKEN_OPAR]:
        rc = False
        break
      node = CNode() 
      node.parent = parent
      if (parent <> None):
        parent.children = parent.children + [node]
      parent = node
      prev_token = TOKEN_OPAR

    elif t == ')':
      if not lookup[prev_token][TOKEN_CPAR]:
        rc = False
        break
      node = node.parent
      parent = node.parent
      prev_token = TOKEN_CPAR

    elif t == ':':
      if not lookup[prev_token][TOKEN_COLON]:
        rc = False
        break
      prev_token = TOKEN_COLON
      pass

    elif t == ';':
      if not lookup[prev_token][TOKEN_SEMICOLON]:
        rc = False
        break
      prev_token = TOKEN_SEMICOLON
      pass

    elif t == ',':
      if not lookup[prev_token][TOKEN_COMMA]:
        rc = False
        break
      prev_token = TOKEN_COMMA
      pass

    else:
      if not lookup[prev_token][TOKEN_STRING]:
        rc = False
        break
      if prev_token == TOKEN_COLON and is_number(t):
        node.length = float(t)
      elif prev_token == TOKEN_CPAR:
        node.label = t
      else:
        node = CNode(t)
        node.parent = parent
        if (parent <> None):
          parent.children = parent.children + [node]
      prev_token = TOKEN_STRING

    if root == None:
      root = node
  
  if prev_token <> TOKEN_SEMICOLON:
    rc = False

  return rc,root
    
def tokenize(s):
  return re.findall("'[^']*'|\"[^\"]*\"|[A-Za-z0-9|+\^\?\._\-]+|\(|\)|,|;|:",s)

def validate_trivial(tokens):
  opar = 0
  cpar = 0
  openpars = 0

  if len(tokens) == 0:
    return False

  if '(' in tokens and tokens[0] <> '(':
    return False

  if (tokens[0] == '('):
    temp = tokens[1:]
    openpars = 1
  else:
    temp = tokens
  
  for t in temp:
    if t == '(':
      if (openpars == 0):
        return False
      opar = opar + 1
      openpars = openpars + 1
    elif t == ')':
      cpar = cpar + 1
      openpars = openpars - 1
    if (openpars < 0):
      return False

  return True

  

def stripwhite(text):
  lstquot = text.split('"')
  for i,p in enumerate(lstquot):
    if not i % 2:
      lstapos = p.split("'")
      for j,q in enumerate(lstapos):
        if not j % 2:
          lstapos[j] = re.sub("\s+", "", q)
      lstquot[i] = "'".join(lstapos)
  return '"'.join(lstquot)
      

def read_file(filename):
  with open(filename) as f:
    content = f.readlines()
  return content

def pynt(filename):
  
  trees = read_file(filename)

  i = 0
  for t in trees:
    
    # TASK 1: Remove whitespace
    #t = stripwhite(t);
    print "\n======= Tree " + str(i) + " ======="
    t = CTree(t)
    t.print_tree()
    print "Leaves: " + str(t.leaves_count)
    print "Inner: " + str(t.inner_count)

    # TASK 2: Parse tree

    #print t.strip();
    #print t
    i = i + 1

def pynt_subtree(filename, list_of_leafs):
  trees = read_file(filename)

  for t in trees:
    t = CTree(t)
    subroot = t.subtree(list_of_leafs)
    return subroot

def pynt_subroot_length(filename, list_of_leafs, key):
  trees = read_file(filename)

  for t in trees:
    t = CTree(t)
    subroot = t.subtree(list_of_leafs)
    print key, subroot.length

def avrg_brln_of_clade(filename, list_of_leafs, key):

  trees = read_file(filename)

  for t in trees:
    t = CTree(t)
    subroot = t.subtree(list_of_leafs)
    leafnodes = []
    tmp_length = 0
    lengths = []

    for leaf in list_of_leafs:
	  leafnodes.append(t.get_node(leaf))

    for node in leafnodes:
      taxon = node.label
      tmp_length = node.length

      while node != subroot:
        node = node.parent
        tmp_length += node.length
      lengths.append(tmp_length)
    print key, float(sum(lengths))/float(len(lengths))

def tip_brln(filename,tip_name):
  trees = read_file(filename)

  for t in trees:
    t = CTree(t)
    brln = t.get_node(tip_name)
    print brln.length

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print " usage: " + sys.argv[0] + " NEWICK-FILE"
    sys.exit(0)

  Porifera = ["Oscarella_carmela","Sycon_ciliatum","Leucosolenia_complicata","Ircinia_fasciculata","Chondrilla_nucula","Ephydatia_muelleri","Amphimedon_queenslandica"]

  avrg_brln_of_clade(sys.argv[1], Porifera, "Porifera")
  pynt_subroot_length(sys.argv[1], Porifera, "Porifera")
