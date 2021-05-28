"""
test_runner_namd.py

test runner_namd.py script(s)
"""

import os
import pytest
import glob
import re

import seekr2.modules.common_base as base
import seekr2.modules.mmvt_base as mmvt_base
import seekr2.modules.common_sim_namd as sim_namd
import seekr2.modules.mmvt_sim_namd as mmvt_sim_namd
import seekr2.modules.runner_namd as runner_namd
import seekr2.tests.make_test_model as make_test_model

def test_Runner_namd_default(tmp_path):
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    mymodel = make_test_model.make_test_model(tmp_path, engine="namd")
    #sim_namd_factory = sim_namd.Sim_namd_factory()
    myanchor = mymodel.anchors[1]
    mmvt_output_filename = os.path.join(
                    tmp_path, myanchor.name, "prod", 
                    "%s%d" % (mmvt_base.NAMDMMVT_BASENAME, 1))
    namd_command = "/home/lvotapka/tmp/namd/Linux-x86_64-g++/namd2"
    namd_arguments = ""
    runner = runner_namd.Runner_namd(mymodel, myanchor, namd_command, 
                                     namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index = \
    runner.prepare()
    #sim_namd_factory = sim_namd.Sim_namd_factory()
    my_sim_namd = mmvt_sim_namd.create_sim_namd(mymodel, myanchor, 
                                                      mmvt_output_filename)
    runner.run(my_sim_namd, mmvt_output_filename+".out")
    myglob = mmvt_output_filename+".out*"
    outfiles = glob.glob(myglob)
    outfile = outfiles[0]
    assert os.path.exists(outfile)
    has_correct_line1 = False
    has_correct_line2 = False
    with open(outfile,"r") as f:
        for line in f:
            print(line)
            if re.search("^WallClock", line):
                has_correct_line1 = True
            if re.search("^WRITING", line):
                has_correct_line2 = True
    assert has_correct_line1 == True
    assert has_correct_line2 == True
    """
    namd_arguments = "+wtf"
    runner = runner_namd.Runner_namd(mymodel, myanchor, namd_command, 
                                     namd_arguments)
    default_output_file, output_basename, state_file_prefix, restart_index = \
    runner.prepare(force_overwrite=True)
    #sim_namd_factory = sim_namd.Sim_namd_factory()
    my_sim_namd = mmvt_sim_namd.create_sim_namd(mymodel, myanchor, 
                                                      mmvt_output_filename)
    runner.run(my_sim_namd, mmvt_output_filename+".out")
    myglob = mmvt_output_filename+".out*"
    outfiles = glob.glob(myglob)
    outfile = outfiles[0]
    assert os.path.exists(outfile)
    has_correct_line1 = False
    has_correct_line2 = False
    with open(outfile,"r") as f:
        for line in f:
            print("line:", line)
            
            if re.search("^WallClock", line):
                has_correct_line1 = True
            if re.search("^WRITING", line):
                has_correct_line2 = True
    
    with pytest.raises(Exception):
        assert has_correct_line1 == True
    with pytest.raises(Exception):
        assert has_correct_line2 == True
    """
    return
