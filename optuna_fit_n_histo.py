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

# get path to folder of this script
current_dir = os.path.dirname(__file__)

# path to theoretical model script
#model_path = "/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/try_stuff/optuna_try/model_try001.C"
model_path = "/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/try_stuff/optuna_try/model_easy_try001.C"
"""
# create study globally
global_directions = ["minimize"]
global_sampler = optuna.samplers.TPESampler(seed=42)
global_study = optuna.create_study(directions=global_directions, sampler=global_sampler) """

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
def worker(cpu, exp_data, par_ranges, path, nbins, trial, chi_sqrs): # used to have study_name and and data_url as input also 
#    study = optuna.load_study(study_name=study_name, storage=storage_url)
    target_function(par_ranges, trial)
    # extract the parameters
    current_parameters = list(trial.params.values())
    # give them to model script
    path_temp_files = path
    executable_path = "/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/CATS/Cproject/bin/LocalFemto"
#    print(current_parameters)
#    print(cpu)
#    print(nbins)
    current_parameters_str = f"{current_parameters}"
    cpu_str = f"{cpu}"
    nbins_str = f"{nbins}"
    command = [executable_path, current_parameters_str, cpu_str, nbins_str]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#    print(f"started process {cpu}...\n")
#    print("communicating...\n")
    stdout, stderr = process.communicate()
    
#    print(stdout)
#    print(stderr)
#    print(process.returncode)
#    print("done communicating :D\n")
    
    # model script returns temporary .root files
    paths = []
    for file in os.listdir(path_temp_files):
#        print("file:", file)
        if file.startswith(f"th_model_cpu{cpu}") and file.endswith(".root"):
            path_j = f"{path}/{file}"
            paths.append(path_j)
#    print("paths:", paths)
    th_histos_names = []
    for i in range(len(list(exp_data.values()))):
        name_i = f"th_histo_{i+1}"
        th_histos_names.append(name_i)
#    print(th_histos_names)
    ordered_paths = []
    for name in th_histos_names:
        for i in range(len(paths)):
            if name in paths[i]:
                ordered_paths.append(paths[i])
    primary_th_dict = dict(zip(th_histos_names, ordered_paths))
#    print("primary dictionary:", primary_th_dict)
    th_histos_dict = get_histos_from_dict(primary_th_dict)
    histos_t = list(th_histos_dict.values())
#    print("histos_t:", histos_t)
    histos_exp = list(exp_data.values())
#    print("histos_exp:", histos_exp)
    # take root files and calculate new chi^2
#    chi_sqrs = []
    for i in range(len(histos_exp)):
        chi_sqr_i = find_chi_sqr(histos_t[i], histos_exp[i])
        chi_sqr_i = chi_sqr_i / nbins[i]
        chi_sqrs.append(chi_sqr_i)
    trial.set_user_attr("cpu", cpu)
#    global_study.tell(trial, chi_sqrs)
    print(f"finished process cpu {cpu}")
    return 0


# new worker after making asking in parent process
def worker_new(cpu, exp_data, path, nbins, trial, param_string, result_queue):
    # give parameters to model script
    path_temp_files = path
    executable_path = "/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/CATS/Cproject/bin/LocalFemto"
    cpu_str = f"{cpu}"
    nbins_str = f"{nbins}"
    command = [executable_path, param_string, cpu_str, nbins_str]
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
    primary_th_dict = dict(zip(th_histos_names, ordered_paths))
    th_histos_dict = get_histos_from_dict(primary_th_dict, all_histos_th_created_bool)
    histos_t = list(th_histos_dict.values())
    histos_exp = list(exp_data.values())
    # take root files and calculate new chi^2
    primary_th_dict = dict(zip(th_histos_names, ordered_paths))
    chi_sqrs = []
    for i in range(len(histos_exp)):
        if histos_t[i] == "no_histo_created":
            print("th histo not made! Setting chi^2 to large value...\n")
            chi_sqrs.append(pow(10, 10))
        else:
            try:
                chi_sqr_i = find_chi_sqr(histos_t[i], histos_exp[i])
                chi_sqr_i = chi_sqr_i / nbins[i]
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

'''
def get_histos_from_dict(histos_dict, all_histos_created_bool):
    histos = {}
    if False not in all_histos_created_bool:
        for name, path in histos_dict.items():
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
    elif False in all_histos_created_bool:
        for name, path in histos_dict.items():
            if path.endswith("_no_such_file.root") == False:
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
'''

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
    if nbins_exp != nbins:
        sys.stderr.write(f"{th_histo} and {exp_histo} have different number of bins! Exiting program...\n")
        sys.exit(1)
    chi_sqr = 0
    for i in range(1, nbins+1):
#        print("th_histo value:", th_histo.GetBinContent(i))
#        print("exp_histo value:", exp_histo.GetBinContent(i))
        value = pow(th_histo.GetBinContent(i) - exp_histo.GetBinContent(i), 2)
#        print("value", value)
        error_i = exp_histo.GetBinError(i)
#        print(error_i)
        chi_sqr += value / pow(error_i, 2)
#        print(chi_sqr)
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
    histogram_names = []
    parameters = []
    th_histo_index = []
    exp_histo_index = []
    cpu_config = 4 # default cpu number :)
#    print(len(config_data[0]))
    for i in range(len(config_data[0])):
        if str(config_data[0][i]).startswith("histo"):
            histogram_names.append(config_data[3][i])
            if config_data[0][i].startswith("histo_th"):
                th_histo_index.append(i)
            if config_data[0][i].startswith("histo_exp"):
                exp_histo_index.append(i)
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
            
#    print(parameters)
    if len(th_histo_index) != len(exp_histo_index):
        sys.stderr.write("Different number theoretical and experimental histograms! Exititng program...\n")
        sys.exit(1)

    par_ranges_df = config_data.iloc[parameters[0]:parameters[-1]+1, 0:3]
    n_data_h = len(histogram_names)
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

# get the theoretical histograms
    th_histo = {} # name of histo - histo
    for i in range(len(th_histo_index)):
        name_i = config_data.iloc[th_histo_index[i], 2]
        path_i = config_data.iloc[th_histo_index[i], 1]
        file_i = ROOT.TFile.Open(path_i, "READ")
        histo = file_i.Get(name_i).Clone()
        if not histo:
            print(f"[ERROR] Could not get histogram '{name_i}' from file '{path_i}'")
            print(f"  - File valid? {file_i.IsZombie()}")
            continue
        histo_i = file_i.Get(name_i).Clone()
        histo_i.SetDirectory(0)
        file_i.Close()
#        print("th name_i", name_i)
 #       print("path_i", path_i)
        name_i = f"teoretical_histo_{i}"
        th_histo[name_i] = histo_i
    th_histo_names_list = list(th_histo.keys())
    print("th histo", th_histo)
    
# get the experimental histograms
    exp_histo = {} # name of histo - histo
    for i in range(len(exp_histo_index)):
        name_i = config_data.iloc[exp_histo_index[i], 2]
        path_i = config_data.iloc[exp_histo_index[i], 1]
        file_i = ROOT.TFile.Open(path_i, "READ")
#        print("exp name_i", name_i)
#        print("path_i", path_i)
        histo_i = file_i.Get(name_i).Clone()
        histo_i.SetDirectory(0)
        file_i.Close()
        name_i = f"experimental_histo_{i}"
        exp_histo[name_i] = histo_i
    exp_histo_names_list = list(exp_histo.keys())

    print("exp histo", exp_histo)

# get chi squared for all histos

#    print(exp_histo)
#    print(th_histo)

    chi_sqr_list = []
    for i in range(len(th_histo_index)):
        th_histo_i = th_histo[th_histo_names_list[i]]
        exp_histo_i = exp_histo[exp_histo_names_list[i]]
        chi_sqr_i = find_chi_sqr(th_histo_i, exp_histo_i)
        print(chi_sqr_i)
        chi_sqr_list.append(chi_sqr_i)
    print("chi_sqr_list:\n", chi_sqr_list)


# create study (abandoned in the name of making it global)
    directions_list_002 = []
    for i in range(len(th_histo_index)):
        directions_list_002.append("minimize")

#    study_name = "Fin01"
    dir_storage_path = os.path.dirname(path_to_config_file)
#    storage_path = f"{dir_storage_path}/optuna_optimizer_data.db"
#    if os.path.exists(storage_path):
#        os.remove(storage_path)
#    data_url = "sqlite:///optuna_optimizer_data.db"
#    data_url = f"sqlite:///{storage_path}"
#    sampler = optuna.samplers.NSGAIISampler()
    sampler = optuna.samplers.TPESampler(seed=1)
#    sampler = optuna.samplers.RandomSampler()
#    study = optuna.create_study(study_name = study_name, storage=data_url, directions=directions_list, sampler=sampler)
    study_002 = optuna.create_study(directions=directions_list_002, sampler=sampler)
    number_of_chi2s = len(exp_histo_names_list)
#    if number_of_chi2s > 1:
#        for i in range (0, number_of_chi2s - 1):
#            global_directions.append("minimize")

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
                p = Process(target=worker_new, args=(cpu, exp_histo, dir_storage_path, nbins_list, trial, param_string, result_queue))
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
            # get the cpu responsible for the best trial in order to delete all other files
            best_cpu = best_trial.user_attrs.get("cpu", None)  # Will be None if not set
            for filename in os.listdir(dir_storage_path):
                old_path = f"{dir_storage_path}/{filename}"
                if filename.startswith(f"th_model_cpu{best_cpu}") and change_best_trial_root_file:
                    print(f"RENAMING cpu {best_cpu} file...\n")
                    filename = "best_trial.root"
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
    
