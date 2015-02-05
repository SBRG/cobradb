# -*- coding: utf-8 -*-

"""
@author: P0N3Y (Pierre SALVY)
Submitted by Amyris (Amyris, Inc). for SBRG@UCSD
"""

from collections import namedtuple
from operator import add, mul, div, sub

'''
Source : http://rosettacode.org/wiki/Parsing/Shunting-yard_algorithm
From here to the end token, this code is a little adaptation of the one mentioned above. They are the ones to
credit for the algorithm.
'''
 
OpInfo = namedtuple('OpInfo', 'prec assoc')
L, R = 'Left Right'.split()
 
ops = {
     '^': OpInfo(prec=4, assoc=R),
     '*': OpInfo(prec=3, assoc=L),
     '/': OpInfo(prec=3, assoc=L),
     '+': OpInfo(prec=2, assoc=L),
     '-': OpInfo(prec=2, assoc=L),
     '(': OpInfo(prec=9, assoc=L),
     ')': OpInfo(prec=0, assoc=L),
     'and': OpInfo(prec=3, assoc=L),
     'or': OpInfo(prec=2, assoc=L),
     }

func_dict = {
     '^': pow,
     '*': mul,
     '/': div,
     '+': add,
     '-': sub,
     #'and': min,
     'or': add,
     }
     
NUM, LPAREN, RPAREN = 'NUMBER ( )'.split()
 
 
def get_input(inp = None):
    'Inputs an expression and returns list of (TOKENTYPE, tokenvalue)'
 
    if inp is None:
        inp = input('expression: ')
    inp = inp.replace('(', ' ( ').replace(')', ' ) ')
    tokens = inp.strip().split()
    tokenvals = []
    for token in tokens:
        if token.lower() in ops:
            tokenvals.append((token, ops[token.lower()]))
        #elif token in (LPAREN, RPAREN):
        #    tokenvals.append((token, token))
        else:    
            tokenvals.append((NUM, token))
    return tokenvals
 
def shunting(tokenvals):
    outq, stack = [], []
    table = ['TOKEN,ACTION,RPN OUTPUT,OP STACK,NOTES'.split(',')]
    for token, val in tokenvals:
        note = action = ''
        if token is NUM:
            action = 'Add number to output'
            outq.append(val)
            table.append( (val, action, ' '.join(outq), ' '.join(s[0] for s in stack), note) )
        elif token.lower() in ops:
            t1, (p1, a1) = token, val
            v = t1
            note = 'Pop ops from stack to output' 
            while stack:
                t2, (p2, a2) = stack[-1]
                if (a1 == L and p1 <= p2) or (a1 == R and p1 < p2):
                    if t1 != RPAREN:
                        if t2 != LPAREN:
                            stack.pop()
                            action = '(Pop op)'
                            outq.append(t2)
                        else:    
                            break
                    else:        
                        if t2 != LPAREN:
                            stack.pop()
                            action = '(Pop op)'
                            outq.append(t2)
                        else:    
                            stack.pop()
                            action = '(Pop & discard "(")'
                            table.append( (v, action, ' '.join(outq), ' '.join(s[0] for s in stack), note) )
                            break
                    table.append( (v, action, ' '.join(outq), ' '.join(s[0] for s in stack), note) )
                    v = note = ''
                else:
                    note = ''
                    break
                note = '' 
            note = '' 
            if t1 != RPAREN:
                stack.append((token, val))
                action = 'Push op token to stack'
            else:
                action = 'Discard ")"'
            table.append( (v, action, ' '.join(outq), ' '.join(s[0] for s in stack), note) )
    note = 'Drain stack to output'
    while stack:
        v = ''
        t2, (p2, a2) = stack[-1]
        action = '(Pop op)'
        stack.pop()
        outq.append(t2)
        table.append( (v, action, ' '.join(outq), ' '.join(s[0] for s in stack), note) )
        v = note = ''
    return table

'''
End of code from http://rosettacode.org/wiki/Parsing/Shunting-yard_algorithm
'''

class Node:
    '''
    Node class to build a tree
    '''
    left , right, data = '', '', 0
    def __init__(self, data):
        # initializes the data members
        self.left = ''
        self.right = ''
        self.data = data
    def __repr__(self):
        return self.data

def parse(expr):
    '''
    Parse from infix notation to postfix(RPN)
    '''
    return shunting(get_input(expr))[-1][2]

def make_tree(rpn):
    '''
    Makes a tree out of a Postfix notation
    '''
    elements = rpn.split()
    stack = []
    for e in elements:
        if e.lower() in ops:
            n = Node(e)
            n.left = stack.pop()
            n.right = stack.pop()
            stack.append(n)
        else:
            stack.append(Node(e))
    return stack

def is_leaf(n):
    '''
    Tells you if it's autumn yet.
    '''
    return True if not n.right and not n.left else False

def getexpr(tree):
    '''
    Transforms a tree into infix notation. Basically your gpr
    '''
    node_value = ""
    if tree:
        node_value = getexpr(tree.left)
        node_value = node_value + ' ' + str(tree.data)
        node_value = node_value + getexpr(tree.right)
        if not is_leaf(tree): node_value = ' (' + node_value + ' )'
    return node_value

def evaluate(tree, and_func = min):
    '''
    Evaluates your tree with the provided AND function. Acceptable AND 
    functions are  min, average, etc ...
    '''
    f = dict(func_dict,**{'and':and_func})
    node_value = 0
    if tree:
        if not is_leaf(tree):
            left = evaluate(tree.left, and_func)
            right = evaluate(tree.right, and_func)
            op = tree.data.lower()
            res = f[op](left,right)
            #print left, op, right, '=',res
            return res
        else:
            return float(tree.data)


if __name__ == '__main__':
    '''
    Honestly, run it once to understand what happens exactly.
    '''
    infix = '((YBR003W) OR (YMR101C) oR (YNR041C) Or (YPL172C) Or (YPR176C AnD YJL031C AND YGL155W))'
    print( 'For infix expression: %r\n' % infix )
    rp = shunting(get_input(infix))
    maxcolwidths = [len(max(x, key=len)) for x in zip(*rp)]
    row = rp[0]
    print( ' '.join('{cell:^{width}}'.format(width=width, cell=cell) for (width, cell) in zip(maxcolwidths, row)))
    for row in rp[1:]:
        print( ' '.join('{cell:<{width}}'.format(width=width, cell=cell) for (width, cell) in zip(maxcolwidths, row)))
    
    print('\n The final output RPN is: %r' % rp[-1][2])

    t = make_tree(rp[-1][2])
    print getexpr(t[0])

    d = {
        'YBR003W':'0',
        'YMR101C':'1',
        'YNR041C':'2',
        'YPL172C':'3',
        'YPR176C':'4',
        'YJL031C':'5',
        'YGL155W':'6',
        }

    s = infix
    for k,v in d.iteritems():
        s = s.replace(k,v)
    print s
    rp = shunting(get_input(s))
    print evaluate(make_tree(rp[-1][2])[0], and_func  = lambda x,y:(x+y)/2)
