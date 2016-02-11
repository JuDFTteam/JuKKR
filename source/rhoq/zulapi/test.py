import subprocess
import os
for i in os.listdir('.'):
    subprocess.call('nm --size-sort --radix=d '+i+'  |  tail -n3',shell=True)

