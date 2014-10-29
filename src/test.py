import re
import subprocess, threading
import os
import signal

class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None
        self.out=u""
        self.err=u""

    def run(self, timeout):
        def target():
            self.process = subprocess.Popen(self.cmd, shell=True,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE, preexec_fn=os.setsid)
            self.out, self.err = self.process.communicate()

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print "killing..."
            self.process.terminate()
            thread.join(1)
            if thread.is_alive():
                print "re-killing..."
                self.process.kill()
                thread.join(1)
                if thread.is_alive():
                    print "re-re-killing..."
                    os.killpg(self.process.pid, signal.SIGTERM) 
            thread.join()

        return (self.process.returncode, self.out, self.err)

cmds = [[u"SIMPLE      SMALL   MPI 1 T 1",
         u"./example/bin/simple -rsa matrix/small.rsa -t 1"],
        [u"SIMPLE      SMALL   MPI 1 T 2",
         u"./example/bin/simple -rsa matrix/small.rsa -t 2"],
        [u"SIMPLE      SMALL   MPI 2 T 1",
         u"mpirun -np 2 ./example/bin/simple -rsa matrix/small.rsa -t 1"],
        [u"SIMPLE      SMALL   MPI 2 T 2",
         u"mpirun -np 2 ./example/bin/simple -rsa matrix/small.rsa -t 2"],
        [u"SIMPLE      ORSIRR  MPI 2 T 2",
         u"mpirun -np 2 ./example/bin/simple -rsa matrix/orsirr.rua -t 2"],
        [u"SIMPLE      YOUNG4C MPI 2 T 2",
         u"mpirun -np 2 ./example/bin/simple -mm matrix/young4c.mtx -t 2"],
        [u"SIMPLE_DIST SMALL   MPI 2 T 2",
         u"mpirun -np 2 ./example/bin/simple_dist -rsa matrix/small.rsa -t 2"],
        [u"SIMPLE_DIST ORSIRR  MPI 2 T 2",
         u"mpirun -np 2 ./example/bin/simple_dist -rsa matrix/orsirr.rua -t 2"],
        [u"SIMPLE_DIST YOUNG4C MPI 2 T 2",
         u"mpirun -np 2 ./example/bin/simple_dist -mm matrix/young4c.mtx -t 2"]

    ]

nb_error = 0
total=0
for cmd in cmds:
    total+=1
    command = Command(cmd[1])
    (ret, out, err) = command.run(timeout=5)
    print cmd[0]+" returned : %d " % ret
    if not ret == 0:
        nb_error+=1
        print "output :"
        for line in out.splitlines():
            print "    ", line
        print "error :"
        for line in err.splitlines():
            print "    ", line

print "%d/%d SUCCESS" %(total-nb_error, total)
if nb_error > 0:
    raise Exception
