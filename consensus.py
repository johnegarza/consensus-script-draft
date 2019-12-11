#This script produces 2-4 files depending on inputs and their contents
#All are packaged together into a folder called hla_calls for convenience
#optitype_calls.txt is always produced, and is essentially a copy of optitype's output
#consensus_calls.txt is also always produced; if no clinical calls are provided, this
#file is identical to optitype_calls.txt. If clinical calls are provided, they are
#reproduced in clinical_calls.txt. If the clinical calls exactly match the optitype calls*, 
#all 3 files described so far will contain the same information, but are not guaranteed to
#be exactly the same (text ordering may differ, depending on the order calls are given in the input). 
#If the clinical calls and optitype calls do not match, mismatched_calls.txt is then produced;
#each line represents a gene. See below (section 'write out call files') for more mismatch details.
#NOTE: optitype only produces MHC class I calls, 

#optitype input format:
#HLA-[A,B,C, ...]*xx:xx
#clinical input format:
# HLA-C*12:02/HLA-C*12:228/HLA-C*12:243 (everything first slash and after is optional)

import sys, os
from collections import defaultdict

####################################
### helper methods for later use ###
####################################

#helper method that takes in the decomposed version of an hla
#string and returns the full delimited string
def build_hla_str(gene, allele_group, spec_allele):
    return gene + "*" + allele_group + ":" + spec_allele

#helper method that takes in a full hla string, like HLA-X*01:02:03:04,
#and splits it into the gene name (HLA-X), allele group (01), and the
#specific allele (02), dropping any fields beyond this, because downstream
#tool do not support these fields
def split_hla_str(full_hla_str):
    gene_name, raw_allele_fields = full_hla_str.split('*')
    split_allele_fields = raw_allele_fields.split(":")
    allele_group_name = split_allele_fields[0]
    specific_allele_name = split_allele_fields[1]
    return (gene_name, allele_group_name, specific_allele_name)

#helper method that creates a mismatch file only if any have been found in the tree,
#and inserts a header upon initially creating the file. Params:
#previously_written- true if the file has already been created; used control header creation
#mismatches- dictionary with sources as keys and a list of alleles called only by that
#            source as values
#returns true if the file was or has ever been written to, false otherwise
def write_mismatch(previously_written, mismatches):
    if (not(mismatches['optitype'] or mismatches['clinical'])):
        #In this case, both arrays are empty, so there's no mismatch to write
        #function has not changed the file state, so return the unmodified flag
        return previously_written

    with open("hla_calls/mismatched_calls.txt", "a") as m_c:
        if not previously_written:
            #add header if this is the first time writing to the file
            m_c.write("optitype_calls\tclinical_calls\n")
        #write the mismatches to the file
        m_c.write( ",".join(mismatches['optitype']) + "\t" + ",".join(mismatches['clinical']) + "\n" )

    return True

########################################
### parse args from the command line ###
########################################

clinical_exists = len(sys.argv) > 2 

optitype_calls = sys.argv[1].split(",")

if clinical_exists:
    raw_clinical_calls = sys.argv[2].split(",")
    #Each clinical call may be a single high confidence call,
    #or a list of uncertain calls separated by slashes
    hc_clinical_calls = []
    u_clinical_calls = []
    for call in raw_clinical_calls:
        if "/" in call:
            u_clinical_calls.append(call)
        else:
            hc_clinical_calls.append(call)

################################################################
### Load HLA types into data structure for consensus calling ###
################################################################

#Create a basic tree out of dictionaries to hold the data from all callers;
#top level keys will be genes, pointing to a nested dictionary with
#allele groups as keys, pointing to a final dictionary with specific alleles
#as keys and the call source as values
# ex: optitype calls HLA-A*01:02 -> {HLA-A: {01: {02: [optitype]}}}

hla_calls = defaultdict( lambda: defaultdict(dict) )

for call in optitype_calls:
    gene, allele_group, spec_allele = split_hla_str(call)

    #checking for existing keys using try/except feels wrong, but follows the python 
    #convention of EAFP (easier to ask forgiveness than permission)
    try:
        #if this path in the tree already exists, a previously-processed identical $call
        #from optitype has already been recorded, so updating the call source is unnecessary
        #this case occurs if an individual is homozygous for a particular gene (according to optitype)
        hla_calls[gene][allele_group][spec_allele]
    except KeyError:
        #if this case is hit, this $call has not been recorded previously and can be safely
        #and normally added to the "tree"
        hla_calls[gene][allele_group][spec_allele] = ['optitype']

if clinical_exists:
    for call in hc_clinical_calls:
        gene, allele_group, spec_allele = split_hla_str(call)

        try:
            #Case 1: this $call was also called by optitype, so add to the 
            #record indicating that clinical data supports this call
            hla_calls[gene][allele_group][spec_allele].append('clinical')
        except KeyError:
            #Case 2: this call is unique to the clinical data, resulting in a
            #KeyError above; create a record and indicate that only clinical data
            #supports this call
            hla_calls[gene][allele_group][spec_allele] = ['clinical']

    for multi_call in u_clinical_calls:
        calls = multi_call.split("/")
        multi_consensus = set()
        for call in calls:
            gene, allele_group, spec_allele = split_hla_str(call)
            try:
                #check if this call already exists in the tree, which will be treated
                #as evidence that this call is the correct call out of the current
                #group of uncertain calls ($multi_call)
                hla_calls[gene][allele_group][spec_allele]

                #this line is only reached if the above has not thrown an exception and
                #moved control to the below except block
                #add as a tuple to avoid re-parsing later
                multi_consensus.add( (gene, allele_group, spec_allele) )
            except KeyError:
                #in this case, there is no further evidence for this particular call
                #nothing to do at this time
                pass

        #if one and only one of the calls from the uncertain group was already in the tree,
        #that is treated as evidence that this particular call was the correct one. It will
        #be accepted and entered into the tree, while the other calls will be discarded
        if len(multi_consensus) == 1:
            accpt_call = multi_consensus.pop()
            hla_calls[accpt_call[0]][accpt_call[1]][accpt_call[2]].append('clinical')
        #otherwise, all uncertain calls from the group will be added to the tree; this means
        #they will be added to the consensus superset (since their validity cannot be disproven),
        #and also used to construct the mismatch file
        else:
            for call in calls:
                gene, allele_group, spec_allele = split_hla_str(call)

                try:
                    hla_calls[gene][allele_group][spec_allele].append('clinical')
                except KeyError:
                    hla_calls[gene][allele_group][spec_allele] = ['clinical']

############################
### write out call files ###
############################

os.mkdir("hla_calls")

#Create an exact copy of optitype calls, to be bundled with other relevant
#files for convenience/later review. Always generated,
with open("hla_calls/optitype_calls.txt", "w") as o_c:
    o_c.write( ",".join(optitype_calls) )

#Create an exact copy of clinical calls, if they exist, to be bundled with 
#other relevant files for convenience/later review.
if clinical_exists:
    with open("hla_calls/clinical_calls.txt", "w") as c_c:
        c_c.write( ",".join(raw_clinical_calls) )

#A consensus file is always generated to be passed on to pvacseq. If there are
#no clinical calls, this file is the same as optitype_calls.txt. If there are, walk
#through the tree and emit everything present as the consensus. If there is a true
#consensus, each class I gene (corresponding to the top level keys of the tree) will have
#at most 2 leaves (1 in the case of a homozygote, or in the rare case that both optitype
#and clinical data only called one allele for this gene), where each leaf represents
#a specific allele call supported by both sources. If there is no true consensus, there
#may be more than 2 leaves per class I gene, and individual leaves may only be supported by
#1 of the 2 sources. These leaves will still be added to the consensus to form a superset,
#since there is not enough evidence to discard them, but they will also be added to a
#mismatch file, which presents side by side lists of the differing alleles called by each
#source, with one gene per line. Note that optitype only makes class I predictions, so any
#class II predictions from the clinical data are always added to the consensus and never
#to the mismatch file
if not clinical_exists:
    with open("hla_calls/consensus_calls.txt", "w") as c_c:
        c_c.write( ",".join(optitype_calls) )
else:
    consensus_calls = []
    mismatch_written = False
    for gene in hla_calls:
        mismatches = {'optitype': [], 'clinical': []}
        for allele_group in hla_calls[gene]:
            for spec_allele in hla_calls[gene][allele_group]:
                callers = hla_calls[gene][allele_group][spec_allele]
                #there are only 3 possibilities for the contents of $callers: it can have
                #one of the two individual caller names, or both names
                if 'clinical' in callers and 'optitype' in callers:
                    consensus_calls.append( build_hla_str(gene, allele_group, spec_allele) )
                else:
                    consensus_calls.append( build_hla_str(gene, allele_group, spec_allele) )
                    mismatches[callers[0]].append( build_hla_str(gene, allele_group, spec_allele) )
        mismatch_written = write_mismatch(mismatch_written, mismatches)

    with open("hla_calls/consensus_calls.txt", "w") as c_c:
        c_c.write( ",".join(consensus_calls) )
