import os,sys
import re
import ROOT
import math
import optuna
import subprocess
import time
import shutil
import multiprocessing
from multiprocessing import Process, Queue
import plotly

import numpy as np
import pandas as pd

from initial_files_needed.AI_scatt_params_to_corr_func import apply_model

# get path to folder of this script
current_dir = os.path.dirname(__file__)

# path to theoretical model script

FmToNu = 5.067731237e-3

start_time = time.time() # start time of the program

# Target function to be optimized by Optuna
# Give values suggested by Optuna to the parameters
def target_function(parameter_ranges, trial):
    sampled_values = {}
    for i in range(len(parameter_ranges)):
        var_name = f"v_par{i+1}"
        var_single_quotes_name = f'v_par{i+1}'
        # Sample the parameters within the specified ranges
        sampled_values[var_name] = trial.suggest_float(var_single_quotes_name, *parameter_ranges[i])
    return 0

# update theoretical histos, calculate new chi^2, update study
def worker_new(cpu, exp_data, path, nbins, trial, param_string, result_queue, k_ranges_list):
    # give parameters to model script
    path_temp_files = path
    executable_path = "/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/CATS/Cproject/bin/LocalFemto"
    cpu_str = f"{cpu}"
    nbins_str = f"{nbins}"
    k_ranges_str = f"{k_ranges_list}"
    command = [executable_path, param_string, cpu_str, nbins_str, k_ranges_str]
#    print("command:\n", command)
    print(f"process cpu {cpu} starting execution of model...")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        print(f"model for cpu {cpu} executed, waiting for stderr...")
        timeout_err = 15 # in seconds
        stdout, stderr = process.communicate(timeout=timeout_err)
        print(f"stderr for cpu {cpu} taken.")
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT={timeout_err}s reached for process cpu {cpu}. Killing process...\n")
        process.kill()
        print(f"killed process cpu {cpu}")
        stdout, stderr = process.communicate()

    print("stdout, stderr:\n", stdout, "\n", stderr)
    th_histos_names = []
    for i in range(len(list(exp_data.values()))):
        name_i = f"th_histo_{i+1}"
        th_histos_names.append(name_i)
        
    # model script returns temporary .root files
    paths = []
    all_histos_th_created_bool = {}
    for file in os.listdir(path_temp_files):
        for name in th_histos_names:
            if file.startswith(f"th_model_cpu{cpu}") and file.endswith(f"{name}.root"):
                path_j = f"{path}/{file}"
                paths.append(path_j)
                all_histos_th_created_bool[name] = True
    if len(th_histos_names) > len(all_histos_th_created_bool):
        for name in th_histos_names:
            if name not in all_histos_th_created_bool.keys():
                all_histos_th_created_bool[name] = False
            
    ordered_paths = []
    for name in th_histos_names:
        for i in range(len(paths)):
            if name in paths[i]:
                ordered_paths.append(paths[i])
#    print("ordered_paths:\n", ordered_paths)
    primary_th_dict = dict(zip(th_histos_names, ordered_paths))
#    print("primary_th_dict:\n", primary_th_dict)
    th_histos_dict = get_histos_from_dict(primary_th_dict, all_histos_th_created_bool)
#    histos_t = list(th_histos_dict.values())
#    histos_exp = list(exp_data.values())
    histos_exp = [exp_data[k] for k in sorted(exp_data)]
    histos_t   = [th_histos_dict[k] for k in sorted(th_histos_dict)]

    # take root files and calculate new chi^2
    primary_th_dict = dict(zip(th_histos_names, ordered_paths))
    chi_sqrs = []

#    print("exp_data:\n", exp_data)
#    print("th_histos_dict:\n", th_histos_dict)
#    print("histos_exp:\n", histos_exp)
#    print("histos_t:\n", histos_t)
#    print("primary th dict:\n", primary_th_dict)
    for i in range(len(histos_exp)):
        if histos_t[i] == "no_histo_created":
            print("th histo not made! Setting chi^2 to large value...\n")
            chi_sqrs.append(pow(10, 10))
        else:
            try:
               # print(f"histos_t[{i}]:", histos_t[i])
               # print(f"histos_exp[{i}]:", histos_exp[i])
                chi_sqr_i = find_chi_sqr(histos_t[i], histos_exp[i])
               # chi_sqr_i = chi_sqr_i / nbins[i]
                chi_sqrs.append(chi_sqr_i)
            except IndexError:
                print("th histo not made! Setting chi^2 to large value...\n")
                chi_sqrs.append(pow(10, 10))
            except:
                print("Some other ERROR occured when calculating CHI^2 in WORKER_NEW function\nExiting...\n")
                return 1
    trial.set_user_attr("cpu", cpu)
    trial_chi_sqrs = [trial, chi_sqrs]
    print(f"putting result cpu {cpu} to queue...\n")
    result_queue.put(trial_chi_sqrs) # return the result via queue
    print(f"result cpu {cpu} has been put in queue.\n")
    print(f"finished process cpu {cpu}")
    
    return 0

# function to get specific histograms from N root files
def get_histos(paths, names_list):
    histos = {}
    for i in range(len(names_list)):
        file_i = ROOT.TFile.Open(paths[i], "READ")
        ROOT.gROOT.cd()
        histo = file_i.Get(names_list[i])
        if not histo:
            print(f"[ERROR] Could not get histogram '{names_list[i]}' from file '{paths[i]}'")
            print(f"  - File valid? {file_i.IsZombie()}")
            continue
        histo_i = file_i.Get(names_list[i]).Clone()
        file_i.Close()
        histos[names_list[i]] = histo_i
    return histos


def get_histos_from_dict(histos_dict, all_histos_created_bool):
    histos = {}
    for name, bool_val in all_histos_created_bool.items():
        if bool_val == True:
            path = histos_dict[name]
            file_i = ROOT.TFile.Open(path, "READ")
            ROOT.gROOT.cd()
            histo = file_i.Get(name)
            if not histo:
                print(f"[ERROR] Could not get histogram '{name}' from file '{path}'")
                print(f"  - File valid? {file_i.IsZombie()}")
                continue
            histo_i = file_i.Get(name).Clone()
            file_i.Close()
            histos[name] = histo_i
        else:
            histos[name] = "no_histo_created"
        
    return histos

def find_chi_sqr(th_histo, exp_histo):
    nbins = th_histo.GetNbinsX()
    nbins_exp = exp_histo.GetNbinsX()
#    print("nbins:", nbins)
#    print("nbins_exp:", nbins_exp)
    if nbins_exp != nbins:
        sys.stderr.write(f"{th_histo} and {exp_histo} have different number of bins! Exiting program...\n")
        sys.exit(1)
    chi_sqr = 0
    min_value = exp_histo.GetXaxis().GetXmin()
    max_value = exp_histo.GetXaxis().GetXmax()
    for i in range(1, nbins+1):
#        print("th_histo value:", th_histo.GetBinContent(i))
#        print("exp_histo value:", exp_histo.GetBinContent(i))
        center_i = exp_histo.GetBinCenter(i)
        if center_i < min_value or center_i > max_value:
            print("HERE'S A BIN OUTSIDE RANGE. SKIPPING...\n")
            continue
        value = pow(th_histo.GetBinContent(i) - exp_histo.GetBinContent(i), 2)
#        print("value", value)
        error_i = exp_histo.GetBinError(i)
#        print("error_i:\n", error_i)
        chi_sqr += value / pow(error_i, 2)
#        print("chi_sqr:\n", chi_sqr)
    return chi_sqr

# check if all chi < chi_target
"""
def check_chi(chi_target, chi_list):
    chi_bool = False
    chi_bool_list = [] # list showing which chis are < target
    bool_list = [] # list with only True to compare
    for i in range(len(chi_list)):
        bool_list.append(True)
        if chi_list[i] < chi_target:
            chi_bool_list.append(True)
        else:
            chi_bool_list.append(False)
    if bool_list == chi_bool_list:
        chi_bool = True
    return chi_bool """


def check_chi(chi_target, chi_list):
    return all(chi < chi_target for chi in chi_list)


def main(path_to_config_file):
# extract paths to all histograms (theoretical data) and their names
    config_data = pd.read_csv(path_to_config_file, sep=r"\s+", header=None)
    config_data.columns = list(range(config_data.shape[1]))
    print(config_data)
    timeout = 10 # default time for calculatino in minutes
# get time limit for calcualtion
    for i in range(len(config_data[0])):
        if config_data[0][i] == "TIMEOUT": 
            try:
                timeout =  float(config_data[1][i])
            except:
                sys.stderr.write("You did not enter a valid number for timeout!\n Using default.\n")
    print("timeout: ", timeout)
# number of histograms to fit and number of parameters
    parameters = []
    th_histo_index = []
    exp_histo_index = []
    k_min_dict = {}
    k_max_dict = {}
    cpu_config = 4 # default cpu number :)
    m_red = 100.0 # reduced mass of the system in MeV
#    print(len(config_data[0]))
    for i in range(len(config_data[0])):
        if config_data[0][i].startswith("histo_exp"):
                exp_histo_index.append(i+1)
                th_histo_index.append(i+1)
        if config_data[0][i].startswith("k_range"):
            try:
                index = int(config_data.iloc[i, 0][-1])
            except:
                sys.stderr.write("Index for k_ranges_index is not integer!\nFix config file.\nExiting program...\n")
                return(1)
            try:
                k_min_dict[index] = float(config_data[1][i])
                k_max_dict[index] = float(config_data[2][i])
            except:
                sys.stderr.write("You did not enter valid k_min and/or k_max!\nFix this.\nExiting program...\n")
                return(1)
        if str(config_data[0][i]).startswith("par"):
#            print(config_data[0][i])
#            print(i)
            parameters.append(i)
        if str(config_data[0][i]) == "CPU":
            try:
                cpu_config = int(config_data[1][i])
            except ValueError:
                print("CPU not integer! Change it to integer in the config file.\n")
                print("Exiting program...\n")
                return(1)
        if str(config_data[0][i]) == "m_red":
            try:
                m_red = float(config_data[1][i])
            except ValueError:
                print("Value for reduced mass is not a number!\nExtiting program...\n")
                sys.exit(1)
            except Error as e:
                print(f"Error: {e}.\n Exiting program...\n")
                sys.exit(1)

    k_ranges_list = []
    keys = list(sorted(k_min_dict.keys()))
    for i in range(len(keys)):
        k_ranges_list.append(k_min_dict[keys[i]])
        k_ranges_list.append(k_max_dict[keys[i]])

    par_ranges_df = config_data.iloc[parameters[0]:parameters[-1]+1, 0:3]
    n_params = len(parameters)
#    print("n_data", n_data_h)
#    print("n_params", n_params)
#    print("parameter ranges\n", par_ranges_df)

# make list of tuples to store the parameter ranges
    par_ranges = []
    for i in range(parameters[0], parameters[-1]+1):
        min_value = float(par_ranges_df[1][i])
        max_value = float(par_ranges_df[2][i])
        tuple_i = (min_value, max_value)
        par_ranges.append(tuple_i)

    
# get the experimental histograms
    exp_histo = {} # name of histo - histo
    th_histo_names_list = []
    for i in range(len(exp_histo_index)):
        name_i = config_data.iloc[exp_histo_index[i]-1, 2]
        path_i = config_data.iloc[exp_histo_index[i]-1, 1]
        file_i = ROOT.TFile.Open(path_i, "READ")
#        print("exp name_i", name_i)
#        print("path_i", path_i)
        histo_i = file_i.Get(name_i).Clone()
        histo_i.SetDirectory(0)
        file_i.Close()
        name_i = f"exp_histo_{i+1}"
        k_min_i = k_min_dict[exp_histo_index[i]]
        k_max_i = k_max_dict[exp_histo_index[i]]
        reduced_bins = histo_i.FindBin(k_max_i) - histo_i.FindBin(k_min_i)
        histo_custom_range_i = ROOT.TH1F(str(name_i), str(name_i), int(reduced_bins), float(k_min_i), float(k_max_i))
#        print(f"k_min_i = {k_min_i}, k_max_i = {k_max_i}\n")
        nbin_min = histo_i.FindBin(k_min_i)
        nbin_max = histo_i.FindBin(k_max_i)
#        print(f"n_min = {nbin_min}, n_max = {nbin_max}\n")
#        print(f"y_kmin = {histo_i.GetBinContent(nbin_min)}\n y_kmax = {histo_i.GetBinContent(nbin_max)}\n")

        j = 0
        for i in range(nbin_min, nbin_max):
            value_i = histo_i.GetBinContent(i)
#            print(f"bin {i}, value = {value_i}\n")
            histo_custom_range_i.SetBinContent(j+1, value_i)
            j = j + 1
        exp_histo[name_i] = histo_custom_range_i
    exp_histo_names_list = list(exp_histo.keys())
    for i in range(len(exp_histo_names_list)):
        name_th_i = f"th_histo_{i+1}"
        th_histo_names_list.append(name_th_i)

    print("exp histo", exp_histo)

# get chi squared for all histos

#    print(exp_histo)
#    print(th_histo)

    chi_sqr_list = []
    for i in range(len(th_histo_index)):
        default_chi_sqr = pow(10, 10)
        chi_sqr_list.append(default_chi_sqr)
    print("chi_sqr_list:\n", chi_sqr_list)


# create study (abandoned in the name of making it global)
    directions_list_002 = []
    for i in range(len(th_histo_index)):
        directions_list_002.append("minimize")

    dir_storage_path = current_dir
    sampler = optuna.samplers.TPESampler(seed=1)
#    sampler = optuna.samplers.RandomSampler(seed=6)
#    sampler = optuna.samplers.NSGAIISampler(seed=9)

    study_002 = optuna.create_study(directions=directions_list_002, sampler=sampler)
    number_of_chi2s = len(exp_histo_names_list)

    current_time = time.time()
    dT = int((current_time - start_time)/60.0)

    chi_target = pow(10, -3)
    chi_bool = check_chi(chi_target, chi_sqr_list)
    print("chi bool: \n", chi_bool)

# get number of bins in each histo :)
    nbins_list = []
    exp_histo_list = list(exp_histo.values())
    for histo in exp_histo_list:
        bins = histo.GetNbinsX()
        nbins_list.append(bins)

    nbins_min = min(nbins_list)

    chi_target = nbins_min
    print(f"chi_target = {chi_target}")

# store best parameters
    best_parameters = []
    best_chi_sqr = list(pow(10, 20) * np.ones(n_params)) # arbitrary big number
    sum_best_chi_sqr = 0
    for i in range(len(best_chi_sqr)):
        sum_best_chi_sqr += best_chi_sqr[i]

    print("timeout:\n", timeout)
    print("dT:\n", dT)
# here - optimize :D
    while dT<timeout:    
        try:
            print("new cycle of while loop...\n")
            if chi_bool == True:
                print("Best parameters:\n", best_parameters)
                print("Best chi^2:\n", best_chi_sqr)
                break
 #           process_trials = {}
            trial_param = {}
            for cpu in range(1, cpu_config+1):
#                trial = global_study.ask()
                trial = study_002.ask()
                trial.set_user_attr("cpu", cpu)
                target_function(par_ranges, trial)
                current_parameters_str = str(list(trial.params.values()))
                current_params_no_m_red = list(trial.params.values())
                current_ss = float(current_params_no_m_red[0])
                current_f = float(current_params_no_m_red[1])
                current_d = float(current_params_no_m_red[2])
                current_a, current_b = apply_model.pass_optuna_trial_guess(m_red, current_f, current_d)
                current_params_list = [m_red, current_ss, float(current_a), float(current_b)]
                current_parameters_str = str(current_params_list)
                print("current_parameters_str:\n", current_parameters_str)
#                process_trials[cpu] = trial
                trial_param_str_list = [trial, current_parameters_str]
                trial_param[cpu] = trial_param_str_list
            print("starting multiple processes...\n")
            result_queue = Queue()
            result_queue.cancel_join_thread()
            processes = []
            chi_sqrs = []
            for cpu, t_p in trial_param.items():
                print(f"starting process: cpu {cpu}")
                trial = trial_param[cpu][0]
                param_string = trial_param[cpu][1]
                p = Process(target=worker_new, args=(cpu, exp_histo, dir_storage_path, nbins_list, trial, param_string, result_queue, k_ranges_list))
                p.start()
               # print("appending process cpu {cpu} in processes...\n")
                processes.append(p)
               # print("process cpu {cpu} appended.\n")

           # for p in processes:
            #    print("do we get here?\n")
             #   p.join() # wait for processes to finish
              #  print("what about here?\n")

            results_expected = cpu_config
            results_collected = 0
            # Collect results from child processes
            while results_collected<results_expected:
                trial, chi_sqrs = result_queue.get()
                study_002.tell(trial, chi_sqrs)
                results_collected +=1

            for i, p in enumerate(processes):
                print(f"Joining process {i}: {p.name}, is_alive: {p.is_alive()}")
                p.join()
                print(f"Joined process {i}")
            

            # Drain queue
            results = []
            while not result_queue.empty():
                try:
                    results.append(result_queue.get_nowait())
                except queue.Empty:
                    break

#            result_queue.close()
 #           print("result_queue is closed\n")
  #          result_queue.join_thread()
   #         print("result_queue joined thread\n")
            
            current_best_params = []
            current_best_chi_sqrs = []
            if study_002.directions and len(study_002.directions) > 1:
                pareto_trials = study_002.best_trials
            else:
                pareto_trials = [study_002.best_trial]
            print("study type:", type(study_002))
            pareto_params = []
            pareto_chi_sqrs = []
            for trial in pareto_trials:
                pareto_params.append(trial.params)
#                print("trial.values:\n", trial.values)
                pareto_chi_sqrs.append(trial.values)
#            print("pareto_chi_sqrs:\n", pareto_chi_sqrs)           
#            chi_2 = [sum(chis_i) for chis_i in pareto_chi_sqrs] # list of chi^2
            chi_2 = []
            for chis_i in pareto_chi_sqrs:
                sum_i = 0
                for i in range(len(chis_i)):
                    sum_i += chis_i[i]
                chi_2.append(sum_i) 
#            print("chi_2:\n", chi_2)
            best_index = chi_2.index(min(chi_2)) # index of smallest sum(chi_i^2)

            current_best_chi_sqrs = pareto_chi_sqrs[best_index]
#            print("pareto trials:", pareto_trials)
            current_best_params = pareto_trials[best_index]
#            print("pareto_trials[best_index]:", current_best_params)
            best_trial = pareto_trials[best_index]

            change_best_trial_root_file = False
            if chi_2[best_index] < sum_best_chi_sqr:
                best_parameters = current_best_params
                best_chi_sqr = current_best_chi_sqrs
                sum_best_chi_sqr = chi_2[best_index]
                change_best_trial_root_file = True
                
            print("current best chi^2: ", chi_2[best_index])
            print("current best parameters:")
            try:
                for key, value in best_trial.params.items():
                        print(f"  {key}: {value}")
            except:
                continue
            # get the cpu responsible for the best trial in order to delete all other files
            best_cpu = best_trial.user_attrs.get("cpu", None)  # Will be None if not set
            for filename in os.listdir(dir_storage_path):
                old_path = f"{dir_storage_path}/{filename}"
                if filename.startswith(f"th_model_cpu{best_cpu}") and change_best_trial_root_file:
                    print(f"RENAMING cpu {best_cpu} file...\n")
                    for th_name in th_histo_names_list:
                        if th_name in filename:
                            filename = f"best_trial_{th_name}.root"
                            new_path = f"{dir_storage_path}/{filename}"
                            os.rename(old_path, new_path)
                if filename.startswith(f"th_model_cpu{cpu}"):
                    print(f"deleting {filename}")
                    os.remove(old_path)
            current_time = time.time()
            dT = int((current_time - start_time)/60.0) # update dT
        except KeyboardInterrupt:
            sys.stderr.write("KeyboardInterrupt! Printing current best parameters...\n")
            break

    print("Best parameters:\n")
    best_params = best_trial.params
    for key, value in best_params.items():
        print(f"  {key}: {value}")
    print("Best chi^2:\n", best_chi_sqr)
    print("chi bool:", chi_bool)
    print("dT:", dT, "min")
    print("deleting unneccesary files...\n")
    for file in os.listdir(dir_storage_path):
        if file.startswith("th_model_cpu"):
            os.remove(f"{dir_storage_path}/{file}")
    print("all junk is deleted\n")
    print("printing last 10 trials to debug:\n")
    for t in study_002.trials[-10:]:
        print(t.number, t.values, t.params)


if __name__ == "__main__":
    path_to_config_file = sys.argv[1]
    multiprocessing.set_start_method('spawn')
    main(path_to_config_file)
    
