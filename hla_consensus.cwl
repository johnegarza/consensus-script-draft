#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: "Script to create consensus from optitype and clinical HLA typing"
requirements:
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "python:3.7.4-slim-buster"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'hla_consensus.py'
        entry: |
                        #This script produces 2-4 files depending on inputs and their contents
                        #All are packaged together into a folder called hla_calls for convenience
                        #optitype_calls.txt is always produced, and is essentially a copy of optitype's output
                        #consensus_calls.txt is also always produced; if no clinical calls are provided, this
                        #file is identical to optitype_calls.txt. If clinical calls are provided, they are
                        #reproduced in clinical_calls.txt. If the clinical calls exactly match the optitype calls*, 
                        #all 3 files described so far will contain the same information, but are not guaranteed to
                        #be exactly the same (text ordering may differ, depending on the order calls are given in the input). 
                        #If the clinical calls and optitype calls do not match*, the clinical calls are treated as the consensus.
                        #mismatched_calls.txt is then produced; each line represents an allele, with the optitype calls on the
                        #left and the differing clinical calls on the right.
                        #*NOTE: optitype only produces MHC class I calls, 

                        import sys, os
    
                        clinical_exists = len(sys.argv) > 2 
                        consensus = []
                        mismatched = False
    
                        optitype_calls = sys.argv[1].split(",")
                        if clinical_exists:
                            clinical_calls = sys.argv[2].split(",")

                            #these will be used to check for mismatches between the clinical and optitype class I calls
                            #optitype results are loaded into *_alleles[0] and clinical data into *_alleles[1]
                            a_alleles = [ set(), set() ]
                            b_alleles = [ set(), set() ]
                            c_alleles = [ set(), set() ]
    
                        os.mkdir("hla_calls")
    
                        with open("hla_calls/optitype_calls.txt", "w") as o_c:
                            o_c.write( ",".join(optitype_calls) )
    
                            for call in optitype_calls:

                                #if clinical data exists, load the optitype data to check for mismatches later
                                #only check for class I alleles because any class II alleles will be from
                                #class II alleles and accepted by default
                                if clinical_exists:
                                    if "HLA-A" in call:
                                        a_alleles[0].add(call)
                                    elif "HLA-B" in call:
                                        b_alleles[0].add(call)
                                    elif "HLA-C" in call:
                                        c_alleles[0].add(call)
    
                        if clinical_exists:
                            with open("hla_calls/clinical_calls.txt", "w") as c_c:
                                c_c.write( ",".join(clinical_calls) )
    
                                for call in clinical_calls:
    
                                    if "HLA-A" in call:
                                        if call in a_alleles[0]:
                                            a_alleles[0].remove(call)
                                        else:
                                            mismatched = True
                                            a_alleles[1].add(call)
                                    elif "HLA-B" in call:
                                        if call in b_alleles[0]:
                                            b_alleles[0].remove(call)
                                        else:
                                            mismatched = True
                                            b_alleles[1].add(call)
                                    elif "HLA-C" in call:
                                        if call in c_alleles[0]:
                                            c_alleles[0].remove(call)
                                        else:
                                            mismatched = True
                                            c_alleles[1].add(call)
    
                                    consensus.append(call)
                        else:
                            consensus = optitype_calls
    
                        with open("hla_calls/consensus_calls.txt", "w") as c_c:
                            c_c.write( ",".join(consensus) )
    
                        if mismatched:
                            with open("hla_calls/mismatched_calls.txt", "w") as m_c:
                                m_c.write( "optitype_calls\tclinical_calls\n" )
                                m_c.write( ",".join(a_alleles[0]) + "\t" + ",".join(a_alleles[1]) + "\n" )
                                m_c.write( ",".join(b_alleles[0]) + "\t" + ",".join(b_alleles[1]) + "\n" )
                                m_c.write( ",".join(c_alleles[0]) + "\t" + ",".join(c_alleles[1]) )            


baseCommand: ['python', 'hla_consensus.py']
inputs:
    optitype_calls:
        type: string[]
        inputBinding:
            position: 1
            itemSeparator: ','
            separate: false
    clinical_calls:
        type: string[]?
        inputBinding:
            position: 2
            itemSeparator: ','
            separate: false
outputs:
    consensus_alleles:
        type: string[]
        outputBinding:
            glob: hla_calls/consensus_calls.txt
            loadContents: true
            outputEval: $(self[0].contents.split(","))
    hla_call_files:
        type: Directory
        outputBinding:
            glob: hla_calls
