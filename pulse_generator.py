import serial
import atexit
import time
import os

__SERIAL = serial.Serial()


def __destructor():
    print('Closing pulse generator connection')

    __close()


atexit.register(__destructor)


def __connect(portno):
    __SERIAL.baudrate = 115200
    __SERIAL.port = 'COM' + str(portno)
    __SERIAL.timeout = 0.03
    __SERIAL.write_timeout = 0.03
    __SERIAL.open()


def __close():
    __SERIAL.close()


def __read_line():
    return __SERIAL.readline()


def connect(portno):
    __connect(portno)
    if __SERIAL.is_open:
        print('pulse generator connected')


def close():
    return __close()


def one_led_pulse(led_num):
    # modified for threading
    # send a pulse in the chosen led number
    if led_num == 0:
        __SERIAL.write("0".encode())
    elif led_num == 1:
        __SERIAL.write("1".encode())
    elif led_num == 2:
        __SERIAL.write("2".encode())
    elif led_num == 3:
        __SERIAL.write("3".encode())
    elif led_num == 4:
        __SERIAL.write("4".encode())
    else:
        __SERIAL.write("0".encode())
    __SERIAL.flush()

def receive():
    # receive the data from the arduino
    line = __SERIAL.readline()
    return line.decode("ascii")
