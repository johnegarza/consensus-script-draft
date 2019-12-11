import subprocess, shutil


optitype_calls = {'HLA-A*01:01', 'HLA-A*01:02'}


def test(clin_calls, expected_consensus):

    #call
    subprocess.call(['python', 'consensus.py', ",".join(optitype_calls), ",".join(clin_calls)])
    
    #get results
    with open('hla_calls/consensus_calls.txt') as f:
        contents = f.read().strip().split(',')
    consensus_result = set(contents)

    #cleanup
    shutil.rmtree('hla_calls')

    #return results
    return expected_consensus == consensus_result


#################
# Run the tests #
#################

clin_and_results = [ ({'HLA-A*01:01', 'HLA-A*01:02'}, {'HLA-A*01:01', 'HLA-A*01:02'}),
                     ({'HLA-A*01:01', 'HLA-A*01:03'}, {'HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03'}),
                     ({'HLA-A*01:03', 'HLA-A*01:04'}, {'HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:04'}),
                     ({'HLA-A*01:01', 'HLA-A*01:01'}, {'HLA-A*01:01', 'HLA-A*01:02'}),
                     ({'HLA-A*01:01', 'HLA-A*02:01'}, {'HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*02:01'}),
                     ({'HLA-A*01:01/HLA-A*01:03', 'HLA-A*01:02'}, {'HLA-A*01:01', 'HLA-A*01:02'}),
                     ({'HLA-A*01:03/HLA-A*01:04', 'HLA-A*01:02'}, {'HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:04'}),
                     ({'HLA-A*01:03/HLA-A*01:04', 'HLA-A*01:05/HLA-A*01:06'}, {'HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:04', 'HLA-A*01:05', 'HLA-A*01:06'}) ]

for input_clin, expected_result in clin_and_results:
    assert( test(input_clin, expected_result) )
