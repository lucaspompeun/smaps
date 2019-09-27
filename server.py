#!/usr/bin/python3
# -*- coding: utf-8 -*-
from flask import Flask, jsonify, send_from_directory, request,render_template, send_file, json
from werkzeug import secure_filename
import os, subprocess, zipfile, tarfile
from main import *
import time
import atexit

from apscheduler.schedulers.background import BackgroundScheduler

app = Flask(__name__, template_folder='template')


class Task(object):
    def __init__(self, task, info):
        self.task = task
        self.info = info
        self.name = info['project']


class Task_Manager(object):
    def __init__(self):
        self.task_list = []
        self.runing_task = []
        self.limit      = 1

    def add_task(self, task, info):
        new_task = Task(task, info)
        self.task_list.append(new_task)
        return True

    def tasks_print(self):
        t = {}
        t['running'] = []
        t['wait_list'] = []
        for task in self.runing_task:
            t['running'].append(task.info)


        for task in self.task_list:
            t['wait_list'].append(task.info)

        return str(json.dumps(t))

    def verif_tasks(self):

        for task in self.runing_task:
            if not task.task.is_alive():
                self.runing_task.remove(task)

        while len(self.runing_task) < self.limit and self.task_list:
            t = self.task_list.pop(0)
            self.runing_task.append(t)
            self.runing_task[-1].task.start()
            print("start task")

        return True


    def num_tasks(self):
        return len(self.task_list) + len(self.runing_task)

    def stop_task(self, task_name):
        print("entrou")
        for task in self.runing_task:
            if task.name == task_name:
                task.task.terminate()

    def move_up(self, name):
        for task in self.task_list:
            if task.name == name:
                i = self.task_list.index(task)
                break

        if i > 0:
            temp = self.task_list[i - 1]
            self.task_list[i -i ] = self.task_list[i]
            self.task_list[i] = temp

    def move_down(self, name):
        for task in self.task_list:
            if task.name == name:
                i = self.task_list.index(task)
                break

        if i < len(self.task_list) -1 :
            temp = self.task_list[i + 1]
            self.task_list[i + i ] = self.task_list[i]
            self.task_list[i] = temp


def get_log(project):
    log_file  = 'projects/' + project + '/log.txt'
    if os.path.isfile(log_file):
        log = open(log_file).readlines()
        log = log[-1].split("\t")[0]
    else:
        log = '0'

    return log

def decommpress(temp_file,orig_file):
    if temp_file.endswith(".gz"):
        tar = tarfile.open(temp_file, "r:gz")
        tar.extractall()
        out = tar.getnames()[0]
        tar.close()
    elif (temp_file.endswith("tar")):
        tar = tarfile.open(temp_file, "r:")
        tar.extractall()
        out = tar.getnames()[0]
        tar.close()

    elif (temp_file.endswith("zip")):
        z = zipfile.ZipFile(temp_file, 'r')
        z.extractall()
        out = z.orig_file
        z.close()

    else:
        out = temp_file


    os.system("mv " + out + " " + orig_file)

    return orig_file

def verify_files(temp_file,orig_file):
    filename = decommpress(temp_file,orig_file)
    os.path.exists(filename)
    return filename

def file_ext(f_type):
    types = {
        'read1' : '.fastq',
        'read2' : '.fastq',
        'reference' : '.fasta',
        #'scaffolds' : '.fasta',
        #'gff' : '.gff'
    }
    return types[f_type]

def get_args_p(xargs):
    x = {'read1' : None ,
        'read2' : None,
        'reference' : None,
        #'scaffolds' : None,
        'project' : None,
        'email' : None,
        #'gff' : None,
        #'gft' : 'contiguator',
        'threads' : '60',
        #'outcontig' : None,
        #'identity' : '100',
        'clusters'  : 5,
        # 'assembler' : 'spades'
        }
    for arg in x:
        if arg in xargs:
            x[arg] = xargs[arg]
    return x

global TASKS
TASKS = Task_Manager()

@app.route('/uploader', methods = ['GET','POST'])
def upload_file():
    xargs = {}
    if request.method == 'POST':
        for f_type in request.files:
            file = request.files[f_type]
            filename = "tmp/" + secure_filename(file.filename)
            file.save(filename)
            filename = verify_files(filename, "tmp/" + request.form['project'] + "_" + f_type + file_ext(f_type))
            xargs[f_type] = filename


        #xargs['assembly'] =  request.form['assembly']
        xargs['project']        = request.form['project']
        xargs['clusters'] = request.form['clusters']
        xargs['email'] =  request.form['email']
        #xargs['gft'] = request.form['gft']
        #xargs['assembler'] = request.form['assembler']

        xargs = get_args_p(xargs)

        print(request.form)
        print(xargs)

        #TASKS.add_task(multiprocessing.Process(target=main,args=(xargs,),),xargs)
        cmd = 'smaps -read1 ' + xargs['read1'] + ' -read2 ' + xargs['read2'] + ' -project ' + xargs['project'] + ' -o ' + xargs['clusters']
        if xargs['reference']:
            cmd += ' -reference  ' + xargs['reference']

        os.system(cmd)

    return 'upload'


@app.route('/status', methods = ['GET','POST'])
def get_status():
    log = get_log(request.form['project'])

    return str(log)


@app.route('/time/<project>', methods = ['GET','POST'])
def get_time(project):

    try:
        x = open('projects/' + project + '/log.txt').readlines()
        x = x[-1].split("\t")[2]
    except:
        x= ['0']


    return x



@app.route('/download/<project>', methods=['GET', 'POST'])
def download(project):
    project_file = 'download/' + project + '.tar.gz'
    print(project_file)
    if os.path.isfile(project_file):

        return send_from_directory(directory='', filename=project_file, as_attachment=True)

    else:
        print('false')
        return 'False'


@app.route('/stop/<project>', methods=['GET', 'POST'])
def kill(project):
    TASKS.stop_task(project)
    return 'kill'

@app.route('/move_down/<project>', methods=['GET', 'POST'])
def move_down(project):
    TASKS.move_down(project)
    return 'down'

@app.route('/move_up/<project>', methods=['GET', 'POST'])
def move_up(project):
    TASKS.move_up(project)
    return 'up'


@app.route('/list', methods=['GET', 'POST'])
def list():
    return TASKS.tasks_print()

@app.route('/server', methods=['GET', 'POST'])
def server_status():
    return "True"

@app.route('/num_tasks', methods=['GET', 'POST'])
def nun_tasks():
    return str(TASKS.num_tasks())

@app.route('/search_project/<project>', methods=['GET', 'POST'])
def search_project(project):
    print('projects' + project)
    return str(os.path.isdir('projects/' + project))

@app.route('/gaps/<project>', methods=['GET', 'POST'])
def gaps(project):
    p = 'projects/' + project + "_a/"
    stats  =  p + "fgap/out.stats"
    gaplog =  p + "combined/gapblaster.log"
    gaplog1 = p + "gapblaster/gapblaster.log"

    Gaps = count_gaps(stats, gaplog, gaplog1)
    return str(Gaps)

@app.route('/config/<project>', methods=['GET', 'POST'])
def get_config(project):

    return send_from_directory(directory='projects/' + project + '/' , filename='conf.json')


    return str(Gaps)

@app.route('/finalized', methods=['GET', 'POST'])
def get_finalized():
    finalized = []
    for p in os.listdir('projects'):
        log = 'projects/' + p + '/log.txt'
        if os.path.isfile(log):
            if int(get_log(p)) == 8:
                conf = 'projects/' + p + '/conf.json'
                if os.path.isfile(conf):
                    finalized.append(json.load(open(conf)))

    return json.dumps(finalized)

scheduler = BackgroundScheduler()
scheduler.add_job(func=TASKS.verif_tasks, trigger="interval", seconds=5)
scheduler.start()
atexit.register(lambda: scheduler.shutdown())

app.run(debug=False,host='0.0.0.0',port='6006')