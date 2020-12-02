from __future__ import print_function
import numpy as np
import os,sys
import cobaya
import pyfisher

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Makes a mock BAO likelihood.')
parser.add_argument("exp_name", type=str,help='Positional arg. Options include boss,boss_data,desi.')
parser.add_argument("--boss-include",     type=str,  default='6df,mgs,lowz,cmass',help="A description.")
parser.add_argument("--input-path",     type=str,  default='input',help="Relative path to directory with experiment info.")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-o","--output",type=str,help="Output root",required=True)
required_args.add_argument("-p","--param-file",type=str,help="Parameter file",required=True)
args = parser.parse_args()

_,fids = pyfisher.get_param_info(args.param_file)

if args.exp_name!='boss_data':
    zs,sig_pers = pyfisher.load_bao_experiment_rs_dV_diagonal(args.exp_name,args.input_path,boss_include=args.boss_include.split(','))
    rs_dv = pyfisher.get_bao_rs_dV(zs,params=fids,engine='camb',de='ppf')
    err = sig_pers * rs_dv / 100.
    cov = np.diag(err**2.)
    np.savetxt(f'{args.output}_{args.exp_name}_mock_bao_covtot.txt',cov,delimiter=' ')
    print(rs_dv,cov)
    with open(f'{args.output}_{args.exp_name}_mock_bao.txt','w') as f:
        for i,z in enumerate(zs):
            f.write(f'{z:.2f} {rs_dv[i]}  rs_over_DV\n')
else:
    zs,ret = pyfisher.get_bao_dr12(params=fids,engine='camb',de='ppf')
    with open(f'{args.output}_{args.exp_name}_mock_bao.txt','w') as f:
        for i,z in enumerate(zs):
            f.write(f'{z:.2f} {ret["DM_over_rs"][i]:.2f}  DM_over_rs\n')
            f.write(f'{z:.2f} {ret["bao_Hz_rs"][i]:.4f}  bao_Hz_rs\n')
    print(ret)

