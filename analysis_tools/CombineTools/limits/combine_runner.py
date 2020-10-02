import subprocess.Popen


class CombineRunner(object):
    def __init__(self, jobs = [], max_jobs=16):
        self.jobs = jobs
        self.max_jobs = max_jobs

    def run_combine(method, b_blind=True):
        """
        spawn parallel combine jobs from datacards
        """

        n_jobs = len(self.jobs)
        n_running = 0
        n_completed = 0
        processes = []
        for signal in signals:        
            b_wait = (n_running >= max_jobs)
            while b_wait:
                n_running = 0
                n_completed = 0
                for proc in processes:
                    if proc[0].poll() == None: n_running += 1
                    else: 
                        n_completed += 1
                        if not proc[1].closed:
                            proc[1].close()
                            print 'Job "{}" has finished.'.format(proc[1].name)
                percent = round(float(n_completed)/float(n_jobs)*100, 1)
                sys.stdout.write( '{0:d} of {1:d} ({2:4.2f}%) jobs done.\r'.format(n_completed, n_jobs, percent))
                sys.stdout.flush()
                time.sleep(10)
                b_wait = (n_running >= max_jobs)

            print 'Spawning job: {}'.format(signal)        
            n_running += 1
            f = open("log/{}.log".format(signal),'w')
            command = "nice -n 10 combine {workspace} -M AsymptoticLimits -n --rMin -2 --rMax 10".format(signal_dir, signal, signal) # --cminDefaultMinimizerTolerance 0
            processes.append((subprocess.Popen(command, stdout=f, shell=True),f))
        b_wait = (n_completed < n_jobs) 
        while b_wait:
            n_running = 0
            n_completed = 0
            for proc in processes:
                if proc[0].poll() == None: n_running += 1
                else: 
                    n_completed += 1
                    if not proc[1].closed:
                        proc[1].close()
                        print 'Job "%s" has finished.' % proc[1].name
            percent = float(n_completed)/float(n_jobs)*100
            sys.stdout.write( '{0:d} of {1:d} ({2:4.2f} %) jobs done.\r'.format(n_completed, n_jobs, percent))
            sys.stdout.flush()
            time.sleep(10)
            b_wait = (n_completed < n_jobs)
        print ''
        print 'Done'

