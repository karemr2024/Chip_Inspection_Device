import sys
from PIL import Image
import PySpin
import time
import myspincam
import pulse_generator as pulser
import os
import scipy.io as sio
import numpy as np
import PWM_Acquisition


def pulse_for_thread(led_num):
    """

    :param led_num: 1: red, 2: green, 3: blue
    :return: return True if led is on False if LED is off.
    """
    if led_num == 0:
        pulser.one_led_pulse(led_num)
        return False
        # return True
    elif led_num == 1:
        pulser.one_led_pulse(led_num)
        return True
    elif led_num == 2:
        pulser.one_led_pulse(led_num)
        return True
    elif led_num == 3:
        pulser.one_led_pulse(led_num)
        return True
    elif led_num == 4:
        pulser.one_led_pulse(led_num)
        return True

def take_img_kill_led(cam: PySpin.CameraPtr, initter, ledon):
    """
    Makes sure LED is on for image to be taken. Takes image as a dictionary containing the
    image_array, image_timestamp, image_ptr. Make sure image is taken before led is turned off.

    :param cam: PySpin camera pointer.
    :param initter: Pre-allocated array according to image size. Use myspincam.initialise_img_arrays(cam).
    :param ledon: True if led is on False if led is off.
    :return: img is an image dictionary containing the image_array, image_timestamp, image_ptr.
    """
    if ledon:
        img, complete = myspincam.get_one_image_array(cam, initter)
        if complete:
            pulse_for_thread(0)
            return img


def pulsing_images_w_pulser(cam: PySpin.CameraPtr, imnum: int):
    """
    This function cycles through the LEDs and takes a picture for each colour.

    :param cam: PySpin camera pointer.
    :param imnum: Number of desired images.
    :return: Dictionary of image lists
    """
    img_dict = dict(red=imnum * [None], green=imnum * [None], blue=imnum * [None]) #, dark=imnum * [None])
    initter = myspincam.initialise_img_arrays(cam)
    cam.BeginAcquisition()
    if cam.EventFrameEnd: #trying out EREN
    # if cam.EventFrameStarts: #trying out EREN
        for i in range(imnum):
            # ledoff = pulse_for_thread(4)
            # img_dict['dark'][i] = take_img_kill_led(cam, initter, ledoff)
            ledon1 = pulse_for_thread(1)
            img_dict['red'][i] = take_img_kill_led(cam, initter, ledon1)
            ledon2 = pulse_for_thread(2)
            img_dict['green'][i] = take_img_kill_led(cam, initter, ledon2)
            ledon3 = pulse_for_thread(3)
            img_dict['blue'][i] = take_img_kill_led(cam, initter, ledon3)
        cam.EndAcquisition()
        return img_dict


def pulse_fast_and_save_pulser(cam, imnum):
    """
    :param cam: PySpin camera pointer.
    :param imnum: Number of desired images.
    :return: PIL images from arrays.

    Each key in the imgs dictionary is a colour. A list of images
    are assigned to each colour, in the order they are taken.
    with one PySPin image dictionary assigned to the keys red, green, blue.
    You may visualise this nested structure as such:

    imgs: dict
    imgs = {'red':[{'image_array':[], 'image_timestamp': int, 'image_ptr': PySpin Image Pointer}, ...],
            'green':[{'image_array':[], 'image_timestamp': int, 'image_ptr': PySpin Image Pointer}, ...],
            'blue':[{'image_array':[], 'image_timestamp': int, 'image_ptr': PySpin Image Pointer}, ...]}
    """
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


def save_as_tiff(r: list, g: list, b: list, imnam: str):
    """
    A folder is created in the desired image name.
    Images are saved as tiff stacks for each colour
    with the number of elements chosen in the experiment.
    This process takes a while.

    :param r: List in which red PIL images are saved.
    :param g: List in which green PIL images are saved.
    :param b: List in which blue PIL images are saved.
    :param imnam: Desired image name.
    """
    print("Saving Images as Tiff Stacks ...")
    t1_start = float(time.clock())
    os.mkdir("C:\Tepegoz\Images" + "\_" + imnam)
    r[0].save("C:\Tepegoz\Images" + "\_" + imnam + "\_" + imnam + "_R.tif",
              compression="tiff_deflate", save_all=True,
              append_images=r[1:])
    g[0].save("C:\Tepegoz\Images" + "\_" + imnam + "\_" + imnam + "_G.tif",
              compression="tiff_deflate", save_all=True,
              append_images=g[1:])
    b[0].save("C:\Tepegoz\Images" + "\_" + imnam + "\_" + imnam + "_B.tif",
              compression="tiff_deflate", save_all=True,
              append_images=b[1:])
    t1_stop = float(time.clock())
    print("pulse start: %f" % t1_start)
    print("pulse stop: %f" % t1_stop)
    print("pulse elapsed: %f" % (t1_stop - t1_start))
    print("Images Successfully Saved as Tiff Stacks.\n")
    pass


def save_as_mat(img_dict: dict, imnam_ipt: str):
    os.mkdir("C:\Tepegoz\Images" + "\_" + imnam_ipt)
    r = {}
    g = {}
    b = {}
    # d = {}
    for i in range(len(img_dict['red'])):
        # d[imnam_ipt + "_" + str(i + 1) + "_D_array"] = np.array(img_dict['dark'][i]['image_array'])
        # d[imnam_ipt + "_" + str(i + 1) + "_D_timestamp"] = img_dict['dark'][i]['timestamp']
        r[imnam_ipt + "_" + str(i + 1) + "_R_array"] = np.array(img_dict['red'][i]['image_array'])
        r[imnam_ipt + "_" + str(i + 1) + "_R_timestamp"] = img_dict['red'][i]['timestamp']
        g[imnam_ipt + "_" + str(i + 1) + "_G_array"] = np.array(img_dict['green'][i]['image_array'])
        g[imnam_ipt + "_" + str(i + 1) + "_G_timestamp"] = img_dict['green'][i]['timestamp']
        b[imnam_ipt + "_" + str(i + 1) + "_B_array"] = np.array(img_dict['blue'][i]['image_array'])
        b[imnam_ipt + "_" + str(i + 1) + "_B_timestamp"] = img_dict['blue'][i]['timestamp']

    # sio.savemat(str("C:\Tepegoz\Images" + "\_" + imnam_ipt + "\D_" + imnam_ipt + ".mat"), d)
    sio.savemat(str("C:\Tepegoz\Images" + "\_" + imnam_ipt + "\R_" + imnam_ipt + ".mat"), r)
    sio.savemat(str("C:\Tepegoz\Images" + "\_" + imnam_ipt + "\G_" + imnam_ipt + ".mat"), g)
    sio.savemat(str("C:\Tepegoz\Images" + "\_" + imnam_ipt + "\B_" + imnam_ipt + ".mat"), b)


def apply_default_experiment_settings(cam):
    """
    Apply default settings. Values are displayed on UI by default.

    :param cam: PySpin camera pointer.
    """
    myspincam.get_cam_and_init(cam)
    myspincam.set_acquisition_mode_continuous(cam)
    # myspincam.set_acquisition_mode_singleframe(cam)
    # myspincam.set_trigger(cam)
    myspincam.set_frame_rate(cam, 1)
    myspincam.disable_auto_exp(cam)
    myspincam.disable_auto_gain(cam)
    myspincam.set_exposure(cam, 4000)
    myspincam.set_gain(cam, 0)
    myspincam.set_gamma(cam, 1)
    myspincam.set_pixel_format(cam)
    myspincam.set_cam_buffer(cam)
    # myspincam.configure_trigger(cam, "SOFTWARE")
    pulser.connect("5")


def apply_custom_experiment_setttings(cam, gain_ipt: float, fps_ipt: int, exposure_ipt: int):
    """
    Apply experiment settings according to user input in UI.
    Currently, trying out different buffer handling modes...

    :param cam: PySpin camera pointer.
    :param gain_ipt:
    :param fps_ipt: Suggested FPS values between 0 & 80. Working on optimisation.
    :param exposure_ipt: Exposure values given in microseconds.
    """
    myspincam.get_cam_and_init(cam)
    myspincam.set_acquisition_mode_continuous(cam)
    myspincam.set_frame_rate(cam, fps_ipt)
    myspincam.disable_auto_exp(cam)
    myspincam.disable_auto_gain(cam)
    myspincam.set_exposure(cam, exposure_ipt)
    myspincam.set_gain(cam, gain_ipt)
    myspincam.set_gamma(cam, 1)
    myspincam.set_pixel_format(cam)
    myspincam.set_cam_buffer(cam)
    # myspincam.configure_trigger(cam, "SOFTWARE")
    pulser.connect("5")


def prepulse(fps_ipt):
    timcon = (1 / fps_ipt)
    prepulse_num = fps_ipt
    for i in range(prepulse_num):
        ledon1 = pulse_for_thread(1)
        time.sleep(timcon)
        ledoff = pulse_for_thread(0)
        ledon2 = pulse_for_thread(2)
        time.sleep(timcon)
        ledoff = pulse_for_thread(0)
        ledon3 = pulse_for_thread(3)
        time.sleep(timcon)
        ledoff = pulse_for_thread(0)


def run_experiment_tiff(imnum_ipt: int, imnam_ipt: str, setting: str, gain_ipt: float, fps_ipt: int, exposure_ipt: int):
    print("Starting experiment")
    __SYSTEM = PySpin.System.GetInstance()
    __CAM = __SYSTEM.GetCameras()[0]

    if setting == "DEFAULT":
        apply_default_experiment_settings(__CAM)
    elif setting == "CUSTOM":
        apply_custom_experiment_setttings(__CAM, gain_ipt, fps_ipt, exposure_ipt)

    time.sleep(3)
    r, g, b = pulse_fast_and_save_pulser(__CAM, imnum_ipt)
    save_as_tiff(r, g, b, imnam_ipt)
    pulser.one_led_pulse(0)
    time.sleep(1)
    __CAM.DeInit()
    pass


def run_experiment_mat(imnum_ipt: int, imnam_ipt: str, setting: str, gain_ipt: float, fps_ipt: int, exposure_ipt: int):
    print("Starting experiment")
    __SYSTEM = PySpin.System.GetInstance()
    __CAM = __SYSTEM.GetCameras()[0]

    if setting == "DEFAULT":
        apply_default_experiment_settings(__CAM)
    elif setting == "CUSTOM":
        apply_custom_experiment_setttings(__CAM, gain_ipt, fps_ipt, exposure_ipt)

    time.sleep(3)
    img_dict = pulsing_images_w_pulser(__CAM, imnum_ipt)
    save_as_mat(img_dict, imnam_ipt)
    pulser.one_led_pulse(0)
    time.sleep(1)
    __CAM.DeInit()
    __SYSTEM.ReleaseInstance()


import os
import shutil

source = 'C:\Tepegoz\Images'
destination = 'C:\Tepegoz\iRiS Kinetics Github'

allfiles = os.listdir(source)

for f in allfiles:
    shutil.move(source + f, destination + f)