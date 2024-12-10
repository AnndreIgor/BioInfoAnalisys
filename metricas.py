import psutil
import time


def get_cpu_model():
    with open('/proc/cpuinfo') as f:
        for line in f:
            if 'model name' in line:
                return line.split(':')[1].strip()

