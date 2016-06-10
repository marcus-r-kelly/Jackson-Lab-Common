import tkinter,tkinter.filedialog
import sys
from lib import markutils as mu
import re
from os.path import isfile
from os import remove
import importlib
from django.forms import model_to_dict

REGEX_COMPILED_TYPE=type(re.compile("foo")) ;
# for type comparisons later on

VALID_FIELDS_INTERACTIONS=["interID", "system", "systemType", "Author", "pmid", "throughput", "score", "modification", "phenotypes", "qualifications", "tags", "srcDB"] 
VALID_FIELDS_INPUTS=["interID", "entrezA", "entrezB", "biogridA", "biogridB", "systematicA", "systematicB", "officialA", "officialB", "synonymsA", "synonymsB", "system", "systemType", "Author", "pmid", "organismA", "organismB", "throughput", "score", "modification", "phenotypes", "qualifications", "tags", "srcDB"] 
VALID_FIELDS_NODES=["entrez","biogrid","systematic","official","synonyms","organism"]

# note that field dictionaries still need to include things like entrezA. HOWEVER, some field information is exclusive to interactions and others are exclusive to nodes

# in different field times, which columns(in order) map to which fields of interaction
# or node data

fdms=( "interID" , "officialA" , "officialB" , "score" , "organismA" , \
    "organismB" , "entrezA" , "entrezB" , "srcDB", "tags" )


fd_biogrid=( "interID" , "entrezA" , "entrezB" , "biogridA" , "biogridB" , \
 "systematicA" , "systematicB" , "officialA" , "officialB" , "synonymsA" , \
 "synonymsB" , "system" , "systemType" , "Author" , "pmid" , "organismA" , \
 "organismB" , "throughput" , "score" , "modification" , "phenotypes" , \
 "qualifications" , "tags" , "srcDB" ) ; 

fd_emili=( "interID" , "entrezA" , "entrezB" , "officialA" , "officialB" ) ; 

fd_chem=("interID","biogridA","entrezA","systematicA","official","synonymsA","organismA","","",\
"Author","pmid","","biogridB","officialB","synonymsB",)
# resume, from biogrid chemical tsv

DEFAULT_FIELD_DICTIONARY =  fdms

keyer=lambda x : x.offical + '_' + x.entrez ;
# function to generate node 'key'
# almost always <symbol>_<entrez id>
# this key is used to look up nodes in node-containing dicts or for
# easy comparison between nodes even if they aren't mapped to the same
# memory object. For example, KRAS_3845 will refer to human KRAS
# in two different datasets even if they are separate dataset objects

FORCE_MATCH_QUAL=True
# I'm not sure why I would ever want this to be false tbh.

debugout=open('.interactors_dbg_log.txt','w') ; 

def ee(ek1,ek2) : 
    """
        Using two edge keys e1 and e2, tests if they are equivalent.
        DIRECTED edges are equivalent only if their keys match.
        UNDIRECTED edges are equivalent if the contained node keys match.
    """

    if ek1 == ek2 : 
        return True ;
    else : 
        nk1_1 = ''
        nk2_1 = ''
        thequal='' ; 

        ek1i=iter(ek1)
        c=next(ek1i) ; 
        while c not in {'>','^'} : 
            nk1_1 += c ; 
            c=next(ek1i) ; 

        edgeicon = c; 
        if edgeicon == '>' : 
            return False ; 
            # means that this is a directed edge, so swapping
            # node key order shouldn't produce a match

        c=next(ek1i)
        # move past edge direction signifier
        while c != ':' : 
            nk2_1 += c ; 
            c=next(ek1i)

        try : 
            c=next(ek1i) ; 
            while True : 
                thequal += c ; 
                c=next(ek1i) ; 
        except StopIteration :
            pass ; 

        if nk2_1 + edgeicon + nk1_1 + ':' + thequal  == ek2 : 
            return True
        else : 
            return False ; 

def ei(ek1) :
    """
        Given an an edge key, return the edge key corresponding
        to an inverted edge. For undirected edges, this returns
        the input.
    """


    # inverts edge if directed, otherwise returns input
    nk1_1 = ''
    nk2_1 = ''
    thequal='' ; 

    ek1i=iter(ek1)
    c=next(ek1i) ; 
    while c not in {'>','^'} : 
        nk1_1 += c ; 
        c=next(ek1i) ; 

    edgeicon = c; 
    if edgeicon == '>' : 
        return ek1 ; 
        # means that this is a directed edge, so swapping
        # node key order shouldn't produce a match

    c=next(ek1i)
    # move past edge direction signifier
    while c != ':' : 
        nk2_1 += c ; 
        c=next(ek1i)

    try : 
        c=next(ek1i) ; 
        while True : 
            thequal += c ; 
            c=next(ek1i) ; 
    except StopIteration :
        pass ; 

    return nk2_1 + edgeicon + nk1_1 + ':' + thequal ; 

class interaction(object):
    """
    Contains data corresponding to observations of protein interactions.
    Fields are basically directly ripped from biogrid. Contains the following
    methods.

    interID---------a unique identifier for the observation
    system----------experimental system
    systemType------physical or genetic (basically not used outside of biogrid)
    author----------author of publication (biogrid only)
    pmid------------pubmed ID of publication
    throughput------high-throughput or low-throughput (biogrid only)
    score-----------in MOST cases, modified NSAF.
    modification----basically ignored
    phenotypes------basically ignored
    qualificiations-in non-biogrid data, used to mark different types of
                    interactions that should NOT be merged into the same
                    edge(e.g. interactions with mutant bait protein)
    tags------------basically ignored
    srcDB-----------source database
    nodeA-----------pointer to node object. For directed nodes, this is the
                    "source." By my convention, in AP/MS the source node is
                    the bait protein.

    nodeB---------- pointer to theother node object.

    """
        
    def __init__(self,\
        # to A from B
        nodeA,
        nodeB,
        directed       = False,
        interID        = "",
        system         = "",
        systemType     = "",
        Author         = "",
        pmid           = "",
        throughput     = "",
        score          = "",
        modification   = "",
        phenotypes     = "",
        qualifications = "",
        tags           = "",
        srcDB          = "" ) :

        self.interID        = interID ; 
        self.system         = system ;
        self.systemType     = systemType ;
        self.Author         = Author ;
        self.pmid           = pmid ;
        self.throughput     = throughput ;
        self.score          = score ;
        self.modification   = modification ;
        self.phenotypes     = phenotypes ;
        self.qualifications = qualifications ;
        self.tags           = tags ;
        self.srcDB          = srcDB ;
        self.nodeA          = nodeA ;
        self.nodeB          = nodeB ;
        self.directed       = directed ; 
        
    def clone(self) :
        """
            Makes a copy of the interaction.
        """
        I=interaction(
        interID        =self.interID,\
        system         =self.system,\
        systemType     =self.systemType,\
        Author         =self.Author,\
        pmid           =self.pmid,\
        throughput     =self.throughput,\
        score          =self.score,\
        modification   =self.modification,\
        phenotypes     =self.phenotypes,\
        qualifications =self.qualifications,\
        tags           =self.tags,\
        srcDB          =self.srcDB,\
        nodeA          =self.nodeA,\
        nodeB          =self.nodeB)

        return I ;

    def __str__( self ):
        return( '; '.join( [
            'interID'        +': '+self.interID,
            'system'         +': '+self.system,
            'systemType'     +': '+self.systemType,
            'Author'         +': '+self.Author,
            'pmid'           +': '+str(self.pmid),
            'throughput'     +': '+self.throughput,
            'score'          +': '+str(self.score),
            'modification'   +': '+self.modification,
            'phenotypes'     +': '+self.phenotypes,
            'qualifications' +': '+self.qualifications,
            'tags'           +': '+self.tags,
            'srcDB'          +': '+self.srcDB,
            'nodeA'          +': '+str(self.nodeA),
            'nodeB'          +': '+str(self.nodeB)] ))
   #def set_nodeA(self,nodeA):
   #    self.nodeA=nodeA ;

   #def set_nodeB(self,nodeB):
   #    self.nodeB=nodeB ;

    def edgekey(self) : 
        """
            Generates key for an edge to which this interaction belongs.
        """
        if self.directed : 
            edgeicon='>' ; 
        else :
            edgeicon='^' ; 
        return self.nodeA.key + edgeicon + self.nodeB.key + ':' + self.qualifications

    def __contains__(self,thing) : 

        if thing in { self.nodeA,self.nodeB,self.nodeA.key,self.nodeB.key } : 
            return True ; 
        else : 
            return False ; 

class bgedge(object):
    """
        Container class of interactions.
        bgedge.to--> destination node for interactions (same conventions as interaction class).
        bgedge.whence--> sourcenode for interactions (i.e. APMS bait)
        bgedge.weight--> # of interactions in edge
        bgedge.interactions--> set of contained interactions
        bgedge.meanscore--> mean score of contained interactions
        bgedge.totalscore--> total score of contained interactions
        bgedge.source--> interaction source string if container edges are uniform,
                         otherwise.
        bgedge.directed--> boolean if edge is directed or not.
        bgedge.qual    --> qualification string of contituent edges

        bgedge.key      --> edge key. Takes the form <nk1>[^>]<nk2>:qual,
                            where the caret or wedge indicates an undirected or directed
                            edge respectively.
    """
    def __init__( self, interaction = None, directed = None, qual = '' ): 

        self.to             = None ; 
        self.whence         = None ; 
        self.weight         = 0 ; 
        self.interactions   = set() ; 
        self.meanscore      = 0.0 ;
        self.totalscore     = 0.0 ;
        self.source         = '' ; 
        self.directed       = None ; 
        self.qual           = qual
        self.key            = '' ; 
        self.p              = 1.0 ; 

        if interaction : 
            self.add_interaction(interaction) ; 

    def add_interaction(self,interaction):
        """
            bgedge.add_interaction(interaction)

            Adds the supplied interaction to this edge's substituent interactions,
            and updates class attributes accordingly.
        """

        if interaction in self.interactions : 
            return ; 

        if ( not self.to and not self.whence ) :
            self.whence         = interaction.nodeA ; 
            self.to             = interaction.nodeB ; 
        elif not { self.to , self.whence } ^ { interaction.nodeA,interaction.nodeB } : 
            # symmetric difference is empty ==> no overlap in nodes on this edge w/ incoming interaction
            pass ; 
        else : 
            raise KeyError('Edge\'s nodes are not empty, but incoming interaction\'s nodes do not match.\n') ; 

        if self.directed is None : 
            self.directed= interaction.directed ; 
        elif self.directed != interaction.directed : 
            raise ValueError('Edge is directed and incoming interaction is not ,or vice versa.\n') ; 

        if self.qual != interaction.qualifications : 
            raise KeyError('Edge\'s and incoming interaction do not have matching qual/qualificaitons.\n') ; 
        

        self.interactions.add(interaction) ; 
        self.weight += 1 ;
        try :
            self.totalscore += float(interaction.score) if interaction.score is not None else 0.0 ; 
            self.meanscore=self.totalscore /self.weight ;
            
        except ValueError :
            pass ;

        if not self.key : 
            if self.directed  : 
                edgeicon='>' 
            else : 
                edgeicon='^'
            self.key = self.whence.key + edgeicon + self.to.key + ':' + self.qual

        self.determine_origins() ;

    def remove_interaction(self,interaction):
        """
            bgedge.remove_interaction(interaction)

            removes substituent interactions from edge,
            and updates attributes accordingly
        """

        if interaction not in self.interactions : 
            return ; 
        self.interactions.remove(interaction) ; 
        self.weight -= 1 ;
        try :
            self.totalscore -= float(interaction.score) ;
            self.meanscore=self.totalscore/self.weight ; 
        except ValueError :
            pass ;

        self.determine_origins() ;


    def nkeys(self) :
        """
            returns a tuple of keys to the nodes
            at either end of this edge
        """
        return (self.whence.key,self.to.key) ;

    def determine_origins(self) :
        """
            changes source attribute depending on substituent
            interaction source strings.
        """

        istr="" ;

        for i in self.interactions :
            if not istr : 
                istr=i.srcDB ;
            elif istr and i.srcDB != istr :
                istr='mixed' ;
                break ; 

        self.source=istr ;

    def connects(self,nodea) :
        """
            If the argument supplied is one of the nodes
            connected by this edge, returns the other node.

            Otherwise, returns None.
        """

        if nodea is self.whence :
            return self.to
        elif nodea is self.to : 
            return self.whence
        else : 
            return None ; 

    def connects_key(self,nodeakey) : 
        """
            If the argument supplied is a key to one of the nodes
            connected by this edge, returns the other node's key.

            Otherwise, returns None.
        """

        if nodeakey == self.whence.key :
            return self.to.key
        elif nodeakey == self.to.key : 
            return self.whence.key
        else : 
            return None ; 

    def sharesends(self,other,directed=False) :
        
        if {other.whence.key,other.to.key} != {self.to.key,self.whence.key} : 
            return False ;
        elif directed and self.to == other.to :
            return True ;
        else : 
            return True ;

   #def sameNode(node1,node2):
   #    if ( node1.official == node2.official and node1.organism == node2.organism ):
   #        return True ;
   #    else: 
   #        return False ;

    def __contains__(self,thing) : 

        if thing in { self.to,self.whence,self.to.key,self.whence.key } : 
            return True ; 
        else : 
            return False ; 

class node(object):
    """
        Contains node (gene) information.
        node.entrez =   entrez ID of node if it corresponds to a real gene ; 
                           otherwise 00
        node.systematic =    biogrid relic
        node.official   =   "official" gene symbol as used by entrez
        node.synonyms   =   "other genes" 
        node.organism   =   taxonomy id (by NCBI conventions, 9606 for humans,
                            10090 for mice

        node.debug      =   flag for debug messages
        node.key        =   key (symbol_entrez id)

        node.interactions   =   set containing interactions that link to this node
        node.edges          =   set containing edges that link to this node
        node.partners       =   set containing nodes that are linked by edges to this
    """

    def __init__( self, official, organism, key, entrez = "", biogrid = "", systematic = "", synonyms = "", debug = False ) :

        self.entrez       = entrez ; 
        self.biogrid      = biogrid ;
        self.systematic   = systematic ; 
        self.official     = official ; 
        self.synonyms     = synonyms ; 
        self.organism     = organism ; #this is the NCBI taxonomy id
        self.debug        = debug ; 
        self.key          = key ; 

        self.interactions = {} ; # key: biogrid id # ; value: interaction
        self.partners     = {} ; # key : partner official name + organism; value : node 
        self.edges        = {} ; # key : partner official name ; value: bgedge

    def __str__( self ):
        return( str(self.key) ) 
    def degree(self,within_edge_set=None) :
        """
            node.degree(self,within_edge_set=None) :

            returns the number of nodes that link to this one.
            If within_edge_set is supplied, ONLY count nodes
            that are linked by one or more of the provided edges

        """
        # python doesn't pass by reference


        if not within_edge_set : 
            return len([ k for k in list(self.edges.keys()) ])
        else : 
            if type(next(iter(within_edge_set))) is str : 
               #return sum([ 1 if k in within_edge_set else 0  for k in self.edges.keys() ]) ;
               return len([ k for k in list(self.edges.keys()) \
                if { e.key for e in self.edges[k] } & within_edge_set  ]) ; 
            else : 
                return sum([ 1 if v & within_edge_set else 0 for v in list(self.edges.values()) ])
            # of nodes (because the edges dict is keyed by node) that have
            # some edges in the acceptable edge set (viz. second line does not produce empty set)

    def binds(self,othernode,within_edge_set=None) : 
        """
            node.degree(self,within_edge_set=None) :

            returns the number of nodes that link to this one.
            If within_edge_set is supplied, ONLY count nodes
            that are linked by one or more of the provided edges

        """

        if type(othernode) is str  : 
            if othernode not in list(self.edges.keys()) : 
                return False ; 
            elif within_edge_set : 
                if type(next(iter(within_edge_set))) is bgedge : 
                    return len( {e for e in self.edges[othernode] } & within_edge_set ) > 0 ;
                else : 
                    return len( { e.key for e in self.edges[othernode] } & within_edge_set) > 0
            else : 
                return True ; 

        elif type(othernode) == type(self) : 
             
            if othernode.key not in list(self.edges.keys()) :
                return False ;
            elif within_edge_set : 
                if type(next(iter(within_edge_set))) is bgedge : 
                    return len( {e for e in self.edges[othernode.key] } & within_edge_set ) > 0 ;
                else : 
                    return len({ e.key for e in self.edges[othernode.key] } & within_edge_set) > 0 ; 
            else : 
                return True ; 
        elif type(othernode) in {list,tuple,set} : 
            return [ self.binds(subnode,within_edge_set=within_edge_set) for subnode in othernode ]
        else : 
            raise TypeError('Argument to binds must be of type str or of type interactors.node\n') ; 


class dataSet(object):

    def __init__( self, nodes = None, default_organism = "9606", i_filter = None, n_filter = None, debug = False,
                  superdebug = False, default_src = "", correction_dict = None, interactions = None, edges = None ):

        self.debug              = debug ; 
        self.infilenames        = list()  ; 
        self.default_organism   = default_organism ; 
        self.i_filter           = i_filter ; # interaction filter? 
        self.n_filter           = n_filter ; # node filter? 
        self.default_src        = default_src
        self.keys               = None ; 
        self.superdebug         = superdebug ;
        if ( self.superdebug ):
            self.debug=True ; 


        # this will again be a set of interactions, but not necessarily with complete data
        if interactions is None :
            self.the_data = dict() ;
        else : 
            self.the_data = interactions ; # allows for a "starting" set of interactions

        if edges is None :
            self.edges    = dict() ;
        else : 
            self.edges    = edges ; # allows for a "starting" set of edges

        if nodes is None :
            self.nodes    = {} ; # empty dictionary
        else : 
            self.nodes    = nodes ; # allows for a "starting" dictionary
            self.keys     = set( self.nodes.keys() ) ; 

        if ( not correction_dict ) :
            self.correction_dict = dict() ;
        else : 
            self.correction_dict = correction_dict ; 

        if ( self.debug ) :
            debugout.write("DEBUG:   Sizes at initialization:\n Interactions : {} ; Nodes : {} ; Corrections : {} \n"\
            .format(len(self.the_data),len(self.nodes),len(self.correction_dict))) ; 


    def _load_data_from_db( self, infobj, logger ):

        from django.core.exceptions import FieldDoesNotExist
        from network.models import Interaction as obj
        from django.db.models import Q
        import operator
        from functools import reduce
        
        infobj      = [ i.upper() for i in infobj ]
        
        # assemble filters
        kwargs_f    = { 'organisma': self.default_organism, 'organismb': self.default_organism }
        kwargs_f[ 'srcdb__in' ] = infobj

        # if we have nodes already, we are looking for interactions between those nodes...
        if self.nodes :
            eids   = [ k.split('_')[1] for k in self.keys ]
            q_list = [Q(x) for x in [ ('entreza__in', eids), ('entrezb__in', eids )]]
            
        kwargs_e    = {}
        alldata     = ''
    
        logger.debug('had q_list')
        if self.i_filter[ 'clude' ] == 'exclude' :
            logger.debug('exclude')
            kwargs_e[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ]
            alldata  = obj.objects.filter( reduce(operator.or_, q_list), **kwargs_f ).exclude( **kwargs_e ).values( 'entreza', 'entrezb' ) 
        elif self.i_filter[ 'clude' ] == 'include' :
            logger.debug('include')
            kwargs_f[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ] 
            alldata  = obj.objects.filter( reduce(operator.or_, q_list), **kwargs_f ).values( 'entreza', 'entrezb' ) 

        logger.debug( len(alldata) )
        #alldata = [ model_to_dict( data ) for data in alldata ]
        logger.debug( len(alldata) )
        eids    = [ [r['entreza'], r['entrezb']] for r in alldata ]
        logger.debug( len(eids) )
        eids    = list(set([ eid for elem in eids for eid in elem ]))
        logger.debug( len(eids) )

        kwargs_f[ 'entreza__in' ] = eids
        kwargs_f[ 'entrezb__in' ] = eids

        # ... but we also need interactions between neighbors
        if self.i_filter[ 'clude' ] == 'exclude' :
            logger.debug('exclude')
            kwargs_e[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ] 
            alldata  = obj.objects.filter( **kwargs_f ).exclude( **kwargs_e )
        elif self.i_filter[ 'clude' ] == 'include' :
            logger.debug('include')
            kwargs_f[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ] 
            alldata  = obj.objects.filter( **kwargs_f )
        
        # we still need to deal with n_filter somehow
        logger.debug( len(alldata) )
        return alldata
        
            
    def load( self, infobj, fd = None, h2m = False, m2h = False, directed = False, qualify = '',
              force_qualify = False, force_score = None, logger = '' ) :
        
        self.infilenames.append( infobj ) ;

        alldata     = self._load_data_from_db( infobj, logger )
        added_iKeys = set(self.the_data.keys()) ;
        added_nKeys = set(self.nodes.keys()) ;
            
        if fd :
            fields_dict = fd ;
        else :
            fields_dict = DEFAULT_FIELD_DICTIONARY ; 
        i = 1
        for data in alldata:
            data   = model_to_dict( data )
            node_was_filtered = False ;
            #logger.debug( data )
            if force_score is not None and data.get('score') is not None : 
                data['score'] = force_score ; 

            #if not i % 100:
            #    logger.debug( i )
            i += 1
            i_keys   = {"a" : "" , "b" : "" }  ; 

            ## nodes must be created before interactions
            for c in ( 'a','b') :

                sym         = data.get( 'symbol' + c )
                entrez      = data.get( 'entrez' + c )
                i_keys[c]   = sym + "_" + str(entrez) ; 

                new_node    = node( entrez       = entrez,
                                    biogrid      = data.get("biogrid"+c),
                                    systematic   = data.get("symbol"+c),
                                    official     = sym,
                                    synonyms     = data.get("synonyms"+c),
                                    # establish use of default organisms where not provided
                                    organism     = data.get( 'organism' + c, self.default_organism ),
                                    key          = i_keys[c] ,
                                    debug        = self.superdebug) ;

                if h2m : 
                    new_node  = mouseify_node(new_node) ; 
                    i_keys[c] = new_node.key ; 

                elif m2h : 
                    new_node  = humanize_node(new_node) ;
                    i_keys[c] = new_node.key ; 

                if ( new_node.key not in added_nKeys ) and ( not self.n_filter or ( self.n_filter and self.n_filter.test(new_node) )):
                    if (self.debug): 
                        debugout.write(\
                         "DEBUG:   Creating node {}, uncreated prior to parsing interaction {}.\n"\
                         .format(new_node.key,data["interid"]))

                    self.nodes.update( { new_node.key : new_node } ) ;
                    added_nKeys.add(new_node.key) ;
                elif ( not self.n_filter or ( self.n_filter and self.n_filter.test(new_node) )) :
                    pass ; 
                else :
                    node_was_filtered = True ;


            if ( node_was_filtered ) :
                continue ; 

            qualify_line = qualify[ data.get( 'srcdb' )]
            if force_qualify : 
                thequal = qualify_line
            else : 
                thequal = data.get( 'qualifications', qualify_line ) ; 
            #logger.debug( 'b: "' + i_keys['b']+'"' )
            #logger.debug( self.nodes )
            if data["interid"] not in added_iKeys :
                thisinteraction = interaction( interID        =    data.get("interid"),\
                                               system         =    data.get("system"),\
                                               systemType     =    data.get("systemtype"),\
                                               Author         =    data.get("author", ''),\
                                               pmid           =    data.get("pmid"),\
                                               throughput     =    data.get("throughput"),\
                                               score          =    data.get("score", ''),\
                                               modification   =    data.get("modification", ''),\
                                               phenotypes     =    data.get("phenotypes", ''),\
                                               qualifications =    thequal,\
                                               tags           =    data.get("tags", ''),\
                                               directed       =    directed,\
                                               srcDB          =    data.get("srcdb"),\
                                               nodeA          =    self.nodes[i_keys['a']],\
                                               nodeB          =    self.nodes[i_keys['b']],) ; 
            else:
                logger.debug("WARNING: Interaction {} on line {} has already been added!\n".format(data["interid"],i))
                continue ;

            self.the_data.update({ thisinteraction.interID : thisinteraction }) ;
            added_iKeys.add(data["interid"]) ;
            #necessary because keys() returns a list, not a set
            
            # if the interaction belongs to a previously uncharacterized edge
            #logger.debug( thisinteraction )
            if thisinteraction.edgekey() not in self.edges and ei(thisinteraction.edgekey()) not in self.edges : 

                # create the edge
                newedge = bgedge(interaction = thisinteraction,\
                                 directed    = thisinteraction.directed,\
                                 qual        = thisinteraction.qualifications)

                # add it to the dataSet's dict of edges
                self.edges.update({ newedge.key : newedge }) ; 

                # connect partners of this edge(may be unecessary) ; 
                newedge.whence.partners.update({ newedge.to.key : newedge.to }) ; 
                newedge.to.partners.update({ newedge.whence.key : newedge.whence }) ; 

                if newedge.whence.key not in newedge.to.edges : 
                    newedge.to.edges.update({ newedge.whence.key : { newedge } }) ; 
                else : 
                    newedge.to.edges[newedge.whence.key].add(newedge) ; 
                    
                if newedge.to.key not in newedge.whence.edges : 
                    newedge.whence.edges.update({ newedge.to.key : { newedge } }) ; 
                else : 
                    newedge.whence.edges[newedge.to.key].add(newedge) ; 

            # if the edge was previously characterized, but backwards
            elif ei(thisinteraction.edgekey()) in self.edges : 
                self.edges[ei(thisinteraction.edgekey())].add_interaction(thisinteraction) ; 
            else :
                self.edges[thisinteraction.edgekey()].add_interaction(thisinteraction) ; 

            thisinteraction.nodeA.interactions.update({ thisinteraction.interID : thisinteraction }) ; 
            thisinteraction.nodeB.interactions.update({ thisinteraction.interID : thisinteraction }) ;

            #self.nodes[thisinteraction.nodeA.key].update_interaction(thisinteraction) ;
            #if ( thisinteraction.nodeA.official != thisinteraction.nodeB.official ):
            #    self.nodes[thisinteraction.nodeB.key].update_interaction(thisinteraction) ;

        logger.debug("Completed.") ;
        
    def parse( self, infobj, sep = '\t', fd = None, h2m = False, m2h = False, directed = False, qualify = '',
               force_qualify = False, force_score = None ) :

        filelines = sum( 1 for line in infobj ) ;
        self.infilenames.append( infobj.name ) ; 
        #sys.stderr.write("DATASET:Opened a file of {} lines.\n".format(filelines));
        sys.stderr.write('Parsing {} ({} lines).\n'.format(infobj.name,filelines)) ; 

        infobj.seek(0) ; # return pointer to beginning of file

        i = 1 ;

        added_iKeys = set(self.the_data.keys()) ;
        added_nKeys = set(self.nodes.keys()) ;

        headerline  = infobj.readline() ;
        headers     = headerline.split(sep) ;

        if fd :
            fields_dict = fd ;
        else :
            fields_dict = DEFAULT_FIELD_DICTIONARY ; 

        dataline = infobj.readline()
        while dataline : 
            # skip all commented lines ; first uncommented line is header
            if dataline[0] != '#' : 
                break ; 

        while dataline : 

            node_was_filtered = False ;
            #if (self.debug):
            #sys.stderr.write("DATASET:Parsed line {: >16} of {: >16}.\r".format(i+1,filelines)) ;
            i += 1;
            mu.waitbar(80 * i / filelines ,80,showPct=True,fill='%')
            sys.stdout.flush() ;

            dataline = dataline.strip() ;
            vals     = dataline.split(sep) ; 
            data     = dict(list(zip(fields_dict,vals))) ; 

            if force_score is not None and data.get('score') is not None : 
                data['score'] = force_score ; 

            orgs     = { "A" : "" , "B" : "" }  ;
            src      = "" 
            i_keys   = {"A" : "" , "B" : "" }  ; 

            # establish use of default organisms where not provided

            ## nodes must be created before interactions
            for c in ( 'A','B') :

                orgs[c]     =   data.get( 'organism' + c, self.default_organism)
                sym         =   data.get( 'official' + c ) ; 
                entrez      =   data.get( 'entrez' + c ) ; 
                i_keys[c]   =   sym + "_" + entrez ; 

                new_node    = node(entrez       = entrez,\
                                   biogrid      = data.get("biogrid"+c),\
                                   systematic   = data.get("systematic"+c),\
                                   official     = sym,\
                                   synonyms     = data.get("synonyms"+c),\
                                   organism     = orgs[c],\
                                   key          = i_keys[c] ,\
                                   debug        = self.superdebug) ;

                if h2m : 
                    new_node  = mouseify_node(new_node) ; 
                    i_keys[c] = new_node.key ; 

                elif m2h : 
                    new_node  = humanize_node(new_node) ;
                    i_keys[c] = new_node.key ; 

                if ( new_node.key not in added_nKeys ) and ( not self.n_filter or ( self.n_filter and self.n_filter.test(new_node) )):
                    if (self.debug): 
                        debugout.write(\
                         "DEBUG:   Creating node {}, uncreated prior to parsing interaction {}.\n"\
                         .format(new_node.key,data["interID"])) ;
                    self.nodes.update( { new_node.key : new_node } ) ;
                    added_nKeys.add(new_node.key) ;
                elif ( not self.n_filter or ( self.n_filter and self.n_filter.test(new_node) )) :
                    pass ; 
                else :
                    dataline=infobj.readline() ;
                    node_was_filtered=True ;


            if ( node_was_filtered ) :
                continue ; 

            src=data.get("srcDB") ; 

            if force_qualify : 
                thequal=qualify
            else : 
                thequal=data.get('qualifications',qualify) ; 

            if data["interID"] not in added_iKeys :
                thisinteraction=interaction(\
                interID        =    data.get("interID"),\
                system         =    data.get("system"),\
                systemType     =    data.get("systemType"),\
                Author         =    data.get("Author"),\
                pmid           =    data.get("pmid"),\
                throughput     =    data.get("throughput"),\
                score          =    data.get("score"),\
                modification   =    data.get("modification"),\
                phenotypes     =    data.get("phenotypes"),\
                qualifications =    thequal,\
                tags           =    data.get("tags"),\
                directed       =    directed,\
                srcDB          =    src,\
                nodeA          =    self.nodes[i_keys['A']],\
                nodeB          =    self.nodes[i_keys['B']],) ; 
            else:
                sys.stderr.write("\rWARNING: Interaction {} on line {} has already been added!\n".format(data["interID"],i)) ;
                dataline=infobj.readline() ;
                continue ;

            if not self.i_filter or ( self.i_filter and self.i_filter.test(thisinteraction) ):
                self.the_data.update({ thisinteraction.interID : thisinteraction }) ;
                added_iKeys.add(data["interID"]) ;
                #necessary because keys() returns a list, not a set
            
                # if the interaction belongs to a previously uncharacterized edge
                if thisinteraction.edgekey() not in self.edges and\
                    ei(thisinteraction.edgekey()) not in self.edges : 

                    # create the edge
                    newedge=bgedge(interaction=thisinteraction,\
                             directed=thisinteraction.directed,\
                             qual=thisinteraction.qualifications)

                    # add it to the dataSet's dict of edges
                    self.edges.update({ newedge.key : newedge }) ; 

                    # connect partners of this edge(may be unecessary) ; 
                    newedge.whence.partners.update({ newedge.to.key : newedge.to }) ; 
                    newedge.to.partners.update({ newedge.whence.key : newedge.whence }) ; 

                    if newedge.whence.key not in newedge.to.edges : 
                        newedge.to.edges.update({ newedge.whence.key : { newedge } }) ; 
                    else : 
                        newedge.to.edges[newedge.whence.key].add(newedge) ; 

                    if newedge.to.key not in newedge.whence.edges : 
                        newedge.whence.edges.update({ newedge.to.key : { newedge } }) ; 
                    else : 
                        newedge.whence.edges[newedge.to.key].add(newedge) ; 

                # if the edge was previously characterized, but backwards
                elif ei(thisinteraction.edgekey()) in self.edges : 
                    self.edges[ei(thisinteraction.edgekey())].add_interaction(thisinteraction) ; 
                else :
                    self.edges[thisinteraction.edgekey()].add_interaction(thisinteraction) ; 

                thisinteraction.nodeA.interactions.update({ thisinteraction.interID : thisinteraction }) ; 
                thisinteraction.nodeB.interactions.update({ thisinteraction.interID : thisinteraction }) ;

                #self.nodes[thisinteraction.nodeA.key].update_interaction(thisinteraction) ;
               #if ( thisinteraction.nodeA.official != thisinteraction.nodeB.official ):
               #    self.nodes[thisinteraction.nodeB.key].update_interaction(thisinteraction) ;
            elif (self.debug):
                debugout.write("DEBUG: Interaction {} excluded after filtering.\n".format(thisinteraction.interID)) ;

            dataline=infobj.readline() ;

        if (self.debug):
            sys.stderr.write("\nCompleted.\n") ;
        else : 
            mu.waitbar(80,80,fill='%',showPct=True) ; 
            sys.stderr.write('\n') ; 
        #infobj.close() ;


    def save(self,f,nodes=None,edges=None):

        if isinstance(f,str) :
            outfile=open(f,'w') ;
        else : 
            outfile=f ;


        if edges is None or type(edges) is not set : 
            edges=set(self.edges.values()) ;
        elif type(next(iter(edges))) is str : 
            edges={ self.edges[ek] for ek in edges if ek in self.edges }
        elif type(next(iter(edges))) == bgedge  : 
            pass ;

        if nodes is None or type(nodes) is not set : 
            nodes=set(self.nodes.values()) ;
        elif type(next(iter(nodes))) is str : 
            nodes={ self.nodes[nk] for nk in nodes if nk in self.nodes }
        elif isinstance(next(iter(nodes)),node) == node  : 
            pass ;



        #f=file(fname,"w") ;
        for n in nodes :
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(n.key,n.entrez,
             n.biogrid,n.systematic,n.official,n.synonyms,n.organism)) ;

        outfile.write("\n") ; # this blank line will signal the difference between nodes and interaction

        interactors_outstring="" ;
        for i in range(0,14) :
            interactors_outstring += "{}\t" ;
        interactors_outstring += "{}\n" ;


        #for i in list(self.the_data.values()) :
        for e in edges : 
            for i in e.interactions :
                if {i.nodeA,i.nodeB}.issubset(nodes) : 
                    dir_icon = '>' if i.directed else '.' ;
                    outfile.write(interactors_outstring.format(i.nodeA.key,i.nodeB.key,\
                     i.interID,i.system,i.systemType,i.Author,i.pmid,i.throughput,\
                     i.score,i.modification,i.phenotypes,i.qualifications,i.tags,\
                     i.srcDB,dir_icon)) ;
                else : 
                    pass ; 

        #f.close() ;

    def load_from(self,infobj) :
        filelines=sum( 1 for line in infobj ) ;
        #sys.stderr.write("DATASET:Opened a file of {} lines.\n".format(filelines));
        infobj.seek(0) ; # return pointer to beginning of file
        i=0 ;

        s=infobj.readline() ;

        # reading nodes
        while s.strip()  :
            node_data=s.strip().split('\t') ;
            if node_data[0] in self.nodes : 
                if self.debug : 
                    sys.stderr.write('DEBUG> node'+node_data[0]+'_'+node_data[1]+' already exists in dataset.') ;
                pass ; 
            else : 
                new_node=node(key=node_data[0],entrez=node_data[1],biogrid=node_data[2],systematic=node_data[3],official=node_data[4],synonyms=node_data[5],organism=node_data[6]) ;

                if not self.n_filter or (self.n_filter and self.n_filter.test(new_node)) :
                    self.nodes.update({ new_node.key : new_node }) ; 

            s=infobj.readline() ;

            i += 1 ; 

            mu.waitbar(80 * i / filelines,80,showPct=True) ; 

        # blank line, False for s.strip() triggers switchover to processing interactions 

        s=infobj.readline() ;
        i += 1 ; 

        added_iKeys=set() ; 
        while s.strip() :
            i_data=s.strip().split('\t') ;

            if not ( self.nodes.get(i_data[0]) and self.nodes.get(i_data[1]) ) : 
                i += 1; 
                s=infobj.readline() ;
                mu.waitbar(80 * i / filelines,80,showPct=True) ; 
                continue ; 

            dir_state=True if i_data[14] == '>' else False ; 

            new_inter=interaction(nodeA=self.nodes[i_data[0]],nodeB=\
             self.nodes[i_data[1]],interID=i_data[2],system=i_data[3],\
             systemType=i_data[4],Author=i_data[5],pmid=i_data[6],throughput=i_data[7],\
             score=i_data[8],modification=i_data[9],phenotypes=i_data[10],\
             qualifications=i_data[11],tags=i_data[12],srcDB=i_data[13],directed=dir_state) ;
        
            # add interaction to info for dataset

            if not self.i_filter or (self.i_filter and self.i_filter.test(new_inter)) : 
                self.the_data.update({ new_inter.interID : new_inter}) ;

                added_iKeys.add(i_data[2]) ;
                #necessary because keys() returns a list, not a set
            
                # if the interaction belongs to a previously uncharacterized edge
                if new_inter.edgekey() not in self.edges and\
                    ei(new_inter.edgekey()) not in self.edges : 

                    # create the edge
                    newedge=bgedge(interaction=new_inter,\
                             directed=new_inter.directed,\
                             qual=new_inter.qualifications)

                    # add it to the dataSet's dict of edges
                    self.edges.update({ newedge.key : newedge }) ; 

                    # connect partners of this edge(may be unecessary) ; 
                    newedge.whence.partners.update({ newedge.to.key : newedge.to }) ; 
                    newedge.to.partners.update({ newedge.whence.key : newedge.whence }) ; 

                    if newedge.whence.key not in newedge.to.edges : 
                        newedge.to.edges.update({ newedge.whence.key : { newedge } }) ; 
                    else : 
                        newedge.to.edges[newedge.whence.key].add(newedge) ; 

                    if newedge.to.key not in newedge.whence.edges : 
                        newedge.whence.edges.update({ newedge.to.key : { newedge } }) ; 
                    else : 
                        newedge.whence.edges[newedge.to.key].add(newedge) ; 

                elif ei(new_inter.edgekey()) in self.edges : 
                    self.edges[ei(new_inter.edgekey())].add_interaction(new_inter) ; 
                else :
                    self.edges[new_inter.edgekey()].add_interaction(new_inter) ; 

                new_inter.nodeA.interactions.update({ new_inter.interID : new_inter }) ; 
                new_inter.nodeB.interactions.update({ new_inter.interID : new_inter }) ;

                #self.nodes[new_inter.nodeA.key].update_interaction(new_inter) ;
               #if ( new_inter.nodeA.official != new_inter.nodeB.official ):
               #    self.nodes[new_inter.nodeB.key].update_interaction(new_inter) ;
            elif (self.debug):
                debugout.write("DEBUG: Interaction {} excluded after filtering.\n".format(new_inter.interID)) ;
            else:
                i += 1; 
                s=infobj.readline() ;
                #mu.waitbar(80 * i / filelines,filelines,showPct=True) ; 
                continue ; 

            # add interaction to info for node A
            # add interaction to info for node B
            #if ( newInter.nodeA != newInter.nodeB )  :
            #    self.nodes[iData[1]].update_interaction(newInter) ;
            
            i += 1; 
            s=infobj.readline() ;
            mu.waitbar(80 * i / filelines,80,showPct=True) ; 

        infobj.close() ;

    def find_offs(self,official) : 

        out=list() ; 

        for n in self.nodes : 
            if n.official == official :
                out.append(n) ; 

        if len(out) == 0 : 
            return None ; 
        elif len(out) == 1 :
            return out.pop()
        else:
            return out ; 

    def clone(self,**kwargs) : 
        from io import StringIO 

        thebuffer=StringIO() ;
        self.save(thebuffer,nodes=kwargs.get('nodes')) ; 
        thebuffer.seek(0)
        other=dataSet(n_filter=kwargs.get('n_filter',self.n_filter),\
                      i_filter=kwargs.get('i_filter',self.i_filter)) ; 
        other.load_from(thebuffer) ; 

        return other ; 

    def __contains__(self,thing) : 

        if thing in self.nodes or thing in self.edges : 
            return True ;
        else : 
            return False ; 


def fetch_nodes(queries,thebiogrid) :

    node_queries=set() ; 

    for query in queries : # these are STRINGS, remember
        
        try : 
            node_queries.add(thebiogrid.nodes[query]) ;
        except KeyError : 
            pass ; 
            #sys.stderr.write("WARNING:    query {} has no node in database.\n".format(query)) ;

    return node_queries ;

def fetch_nodes_list(queries,thebiogrid) : 

    node_queries=list() ;

    for query in queries : 
        node_queries.append(thebiogrid.nodes.get(query,None)) ;

    return node_queries ; 
    

def fetch_nodes_regex(queries,thebiogrid,field="official") :

    qtype=type(next(iter(queries))) ; 
    regexset=set() ;

    if qtype is str :
        for query in queries : 
            regexset.add(re.compile(query,re.IGNORECASE)) ; 

    elif qytpe is REGEX_COMPILED_TYPE  :
        for query in queries : 
            regexset.add(query) ;
    else:
        sys.stderr.write("Container for queries does not hold items of type str or compiled regex.\n") ;
        raise ValueError

    node_queries=set() ;

    search_set=set(thebiogrid.nodes.values()) ;

    while search_set :
        s=search_set.pop() ;

        for r in regexset : 
            if r.match(s.__getattribute__(field)) :
                node_queries.add(s) ;
                break ;

    return node_queries ; 
            
def print_edges( edges,print_headers=True,print_mean_scores=False,print_weights=False,print_total_scores=False,\
 sep='\t',fname="",print_organisms=False,print_source=False,inter_string='pp',transform_scores=None,\
 print_quals=True,print_pps=False):

    from numpy import log10

    if ( not fname ) :
        f=sys.stdout  ;
    #elif isfile(fname): 
        #sys.stdout.write('Appending to existing file {}\n'.format(fname)) ; 
        #f=open(fname,"a") ; 
    else : 
        f=open(fname,"w") ; 

    if ( print_headers) :
        f.write("Left{}Inter{}Right".format(sep,sep)) ;
        if ( print_quals ) :
            f.write("{}Qual".format(sep)) ;
        if ( print_weights ) :
            f.write("{}Wt".format(sep)) ;
        if ( print_mean_scores ) :
            f.write("{}Mean".format(sep)) ;
        if ( print_total_scores) :
            f.write("{}Total".format(sep)) ;
        if ( print_organisms ) :
            f.write("{}OrgA{}OrgB".format(sep,sep)) ; 
        if ( print_source ) :
            f.write("{}Source".format(sep,sep)) ; 
        if print_pps : 
            f.write("{}logp").format(sep) ; 

        f.write("\n") ;


    for edge in edges : 

        #edget=tuple(edge.nodes) ;
        if edge.directed : 
            thedir = '>' ; 
        else : 
            thedir = '^' ; 

        f.write("{}{}{}{}{}{}".format(edge.whence.official,sep,thedir,inter_string,sep,edge.to.official)) ;

        if ( print_quals ):
            f.write("{}{}".format(sep,edge.qual)) ;
        if ( print_weights ):
            f.write("{}{}".format(sep,edge.weight)) ;
        if ( print_mean_scores ):
            if transform_scores : 
                f.write("{}{}".format(sep,transform_scores(edge.meanscore))) ;
            else : 
                f.write("{}{}".format(sep,edge.meanscore)) ;
        if ( print_total_scores):
            if transform_scores : 
                f.write("{}{}".format(sep,transform_scores(edge.totalscore))) ;
            else : 
                f.write("{}{}".format(sep,edge.totalscore)) ;
        if ( print_organisms):
            if ( len(edge.nodes) == 2 ):
                f.write("{}{}{}{}".format(sep,edget[0].organism,sep,edget[1].organism)) ;
            else :
                f.write("{}{}{}{}".format(sep,edget[0].organism,sep,edget[0].organism)) ;
        if ( print_source):
            f.write("{}{}".format(sep,edge.source)) ;
        if print_pps : 
            f.write("{}{}".format(sep,-1*log10(edge.p)))

        f.write("\n") ;

    f.close()



def humanize_node(n) :

    if n.organism != '10090' : 
        return n ; 

    from lib import rbase
    rbase.load('m2h') ; 
    rbase.load('hmg') ;

    if rbase.m2h.get(n.entrez) is None :
        h_entrez= n.entrez  ; 
    else :
        try : 
            for e in rbase.m2h[n.entrez] :
                if rbase.hmg['EID'][e]['Symbol'] == n.official.upper() : 
                    h_entrez=e ; 
                    break ; 
            else : 
                h_entrez=next(iter(rbase.m2h[n.entrez])) ; 

            if len(rbase.m2h[n.entrez]) > 1  :
                sys.stderr.write('\rWARNING: Interactors.humanize_node :'\
                +' {}:{} is not the unique homolog of {}:{}.\n'\
                .format(h_entrez,rbase.hmg['EID'][h_entrez]['Symbol'],\
                n.entrez,n.official))
        except KeyError : 
            return n ; 


    try : 
        n.official=rbase.hmg['EID'][h_entrez]['Symbol'] ; 
        n.organism='9606' ;
        n.entrez=h_entrez
        n.key=n.official+'_'+n.entrez ; 
    except(KeyError) : 
        return n ; 

    return n ; 
