import sys
from PIL import Image
import PySpin
import time
import myspincam
import pulse_generator as pulser
import os

# import concurrent.futures
# from threading import Event
# import numpy as np
# import itertools
# import threading
# import cv2
# import serial
# import functools


sys.setrecursionlimit(1000)

__SYSTEM = PySpin.System.GetInstance()
__CAM = __SYSTEM.GetCameras()[0]


def pulse_for_thread(led_num):
    print("ledon")
    if led_num == 0:
        pulser.one_led_pulse(led_num)
        print("0")
    elif led_num == 1:
        pulser.one_led_pulse(led_num)
        print("1")
    elif led_num == 2:
        pulser.one_led_pulse(led_num)
        print("2")
    elif led_num == 3:
        pulser.one_led_pulse(led_num)
        print("3")
    return True


def take_img_kill_led(cam: PySpin.CameraPtr, initter, ledon):
    if ledon:
        img, complete = myspincam.get_one_image_array(cam, initter)
        print("img")
        if complete:
            pulse_for_thread(0)
    return img


def pulsing_images_w_pulser(cam: PySpin.CameraPtr, imnum: int):
    bum = dict(red=imnum * [None], green=imnum * [None], blue=imnum * [None])
    initter = myspincam.initialise_img_arrays(cam)
    for i in range(imnum):
        ledon1 = pulse_for_thread(1)
        bum['red'][i] = take_img_kill_led(cam, initter, ledon1)
        ledon2 = pulse_for_thread(2)
        bum['green'][i] = take_img_kill_led(cam, initter, ledon2)
        ledon3 = pulse_for_thread(3)
        bum['blue'][i] = take_img_kill_led(cam, initter, ledon3)
    return bum


def pulse_fast_and_save_pulser(cam, imnum):
    t1_start = float(time.clock())
    imgs = pulsing_images_w_pulser(cam, imnum)
    t1_stop = float(time.clock())
    print("pulse start: %f" % t1_start)
    print("pulse stop: %f" % t1_stop)
    print("pulse elapsed: %f" % (t1_stop - t1_start))
    print("Saving Images as Arrays...")
    r = []
    g = []
    b = []
    for i in range(len(imgs['red'])):
        # imgs[i]['red']['image_ptr'].Save(str(i) + "_R_" + str(imnam) + ".tif")
        # Keep in mind that imnam has been erased as input argument.
        r.append(Image.fromarray(imgs['red'][i]['image_array']))
        # print(imgs['red'][i]['timestamp'])

        # imgs[i]['green']['image_ptr'].Save(str(i) + "_G_" + str(imnam) + ".tif")
        # Keep in mind that imnam has been erased as input argument.
        g.append(Image.fromarray(imgs['green'][i]['image_array']))
        # print(imgs['green'][i]['timestamp'])

        # imgs[i]['blue']['image_ptr'].Save(str(i) + "_B_" + str(imnam) + ".tif")
        # Keep in mind that imnam has been erased as input argument.
        b.append(Image.fromarray(imgs['blue'][i]['image_array']))
        # print(imgs['blue'][i]['timestamp'])

    return r, g, b


def save_as_tiff(r: list, g: list, b:list, imnam: str):
    print("Saving Images as Tiff Stacks ...")
    t1_start = float(time.clock())
    os.mkdir("C:\Tepegoz\Images" + "\_" + imnam)
    r[0].save("C:\Tepegoz\Images" + "\_" + imnam + "_R.tif", compression="tiff_deflate", save_all=True, append_images=r[1:])
    g[0].save("C:\Tepegoz\Images" + "\_" + imnam + "_G.tif", compression="tiff_deflate", save_all=True, append_images=g[1:])
    b[0].save("C:\Tepegoz\Images" + "\_" + imnam + "_B.tif", compression="tiff_deflate", save_all=True, append_images=b[1:])
    t1_stop = float(time.clock())
    print("pulse start: %f" % t1_start)
    print("pulse stop: %f" % t1_stop)
    print("pulse elapsed: %f" % (t1_stop - t1_start))
    print("Images Successfully Saved as Tiff Stacks.\n")
    pass


def apply_default_experiment_settings(cam):
    myspincam.get_cam_and_init(cam)
    myspincam.set_acquisition_mode_continuous(cam)
    # myspincam.set_acquisition_mode_singleframe(cam)
    # myspincam.set_trigger(cam)
    myspincam.disable_auto_exp(cam)
    myspincam.disable_auto_gain(cam)
    myspincam.set_exposure(cam, 4000)
    myspincam.set_gain(cam, 0)
    myspincam.set_gamma(cam, 1)
    myspincam.set_pixel_format(cam)
    pulser.connect("5")


apply_default_experiment_settings(__CAM)

imnam_ipt = input("Input experiment details in format Xnm_Xms_MMDDYY_Camera: ")
imnum_ipt = int(input("Number of pulses: "))
prepulse_ipt = int(input("Number of prepulses: "))
choice_ipt = input("Type 'E' for experimental settings, otherwise type 'N'.")

__CAM.BeginAcquisition()
# pulsing_LED_fast(__ARDUINO, prepulse_ipt)
t1_start = float(time.clock())
r, g, b = pulse_fast_and_save_pulser(__CAM, imnum_ipt)
save_as_tiff(r, g, b, imnam_ipt)
pulser.one_led_pulse(0)
t1_stop = float(time.clock())
time.sleep(1)
__CAM.EndAcquisition()
__CAM.DeInit()

print("pulse start: %f" % t1_start)
print("pulse stop: %f" % t1_stop)
print("pulse elapsed: %f" % (t1_stop - t1_start))
