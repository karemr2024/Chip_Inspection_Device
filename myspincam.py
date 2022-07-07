import PySpin
import numpy as np


def get_cam_and_init(cam):
    if cam is None:
        raise RuntimeError('Camera not found')
    else:
        assert isinstance(cam, PySpin.CameraPtr)
        cam.Init()
        return cam


def initialise_img_arrays(cam):
    nodemap = cam.GetNodeMap()
    node_height = PySpin.CIntegerPtr(nodemap.GetNode('Height'))
    node_width = PySpin.CIntegerPtr(nodemap.GetNode('Width'))
    image_size = [node_height.GetValue(), node_width.GetValue()]
    imray_initialiser = np.ones((image_size[0], image_size[1]), dtype=np.uint16)
    return imray_initialiser


def get_one_image_fast(cam):
    imdict = {}
    img = cam.GetNextImage()
    while img.IsIncomplete():
        img = cam.GetNextImage()
    if not img.IsIncomplete():
        # imray = np.array(img.GetNDArray(), dtype=np.uint16)
        imray = np.frombuffer(img.GetNDArray(), dtype=np.uint16)
        imdict['image_ptr'] = img
        imdict['image_array'] = imray
        imdict['timestamp'] = img.GetTimeStamp()
        img.Release()
        return imdict


def get_one_image_array(cam, imray):
    imdict = {}
    if cam.EventFrameStart: #original
    # if cam.EventFrameEnd:
    # if cam.EventExposureEnd:
    # if cam.EventExposureStart:
        img = cam.GetNextImage()
        while img.IsIncomplete():
            img = cam.GetNextImage()
            img.Release()
        if not img.IsIncomplete():
            # if cam.EventExposureEnd: #original
            # if cam.EventFrameStart:
            if cam.EventFrameEnd:
            # if cam.EventExposureStart:
                img.Release()
                imray = np.multiply(imray, np.array(img.GetNDArray(), dtype=np.uint16))
                # imray = np.multiply(imray, np.frombuffer(img.GetNDArray(), dtype=np.uint16))
                imdict['image_ptr'] = img
                imdict['image_array'] = imray
                imdict['timestamp'] = img.GetTimeStamp()
                return imdict, True


def cam_node_cmd(cam, cam_attr_str, cam_method_str, pyspin_mode_str=None, cam_method_arg=None):
    # Performs cam_method on input cam with optional access mode check
    # First, get camera attribute
    cam_attr = cam
    cam_attr_str_split = cam_attr_str.split('.')
    for sub_cam_attr_str in cam_attr_str_split:
        cam_attr = getattr(cam_attr, sub_cam_attr_str)

    # Print command info
    info_str = 'Executing: "' + '.'.join([cam_attr_str, cam_method_str]) + '('
    if cam_method_arg is not None:
        info_str += str(cam_method_arg)
    print(info_str + ')"')

    # Perform optional access mode check
    if pyspin_mode_str is not None:
        if cam_attr.GetAccessMode() != getattr(PySpin, pyspin_mode_str):
            raise RuntimeError('Access mode check failed for: "' + cam_attr_str + '" with mode: "' +
                               pyspin_mode_str + '".')

    # Format command argument in case it's a string containing a PySpin attribute
    if isinstance(cam_method_arg, str):
        cam_method_arg_split = cam_method_arg.split('.')
        if cam_method_arg_split[0] == 'PySpin':
            if len(cam_method_arg_split) == 2:
                cam_method_arg = getattr(PySpin, cam_method_arg_split[1])
            else:
                raise RuntimeError('Arguments containing nested PySpin arguments are currently not '
                                   'supported...')

    # Perform command
    if cam_method_arg is None:  # pylint: disable=no-else-return
        return getattr(cam_attr, cam_method_str)()
    else:
        return getattr(cam_attr, cam_method_str)(cam_method_arg)


def cam_node_editor(cam, cam_attr_str, __cam_node_cmd, cam_method_str, pyspin_mode_str=None, cam_method_arg=None):
    return cam_node_cmd(cam, get_cam_and_init(cam), cam_attr_str, cam_method_str, pyspin_mode_str, cam_method_arg)


def set_acquisition_mode_continuous(cam):
    cam_node_cmd(cam, 'AcquisitionMode', 'SetValue', 'RW', PySpin.AcquisitionMode_Continuous)
    print('Setting Acquisition Mode to Continuous')


def set_acquisition_mode_singleframe(cam):
    cam_node_cmd(cam, 'AcquisitionMode', 'SetValue', 'RW', PySpin.AcquisitionMode_SingleFrame)
    print('Setting Acquisition Mode to Single Frame')


def set_gain(cam, gain):
    cam_node_cmd(cam, 'Gain', 'SetValue', 'RW', gain)
    print('Setting Gain to %d' % gain)


def set_exposure(cam, exposure):
    cam_node_cmd(cam, 'ExposureTime', 'SetValue', 'RW', exposure)
    print('Setting Exposure to %d' % exposure)


def set_frame_rate(cam, frame_rate):
    cam_node_cmd(cam, 'AcquisitionFrameRateEnable', 'SetValue', 'RW', True)
    cam_node_cmd(cam, 'AcquisitionFrameRate', 'SetValue', 'RW', frame_rate)


def set_gamma(cam, gamma_val):
    print('Setting Gamma to ' + str(gamma_val))
    cam_node_cmd(cam, 'Gamma', 'SetValue', 'RW', gamma_val)


def set_pixel_format(cam):
    # Sets pixel format of the camera
    if cam.PixelFormat.GetAccessMode() != PySpin.RW:
        print('Unable to set pixel format. Aborting...')
        return False

    cam.PixelFormat.SetValue(PySpin.PixelFormat_Mono16)
    print('Pixel format is set %s...' % cam.PixelFormat.GetCurrentEntry().GetSymbolic())


def set_trigger(cam):
    cam_node_cmd(cam, 'TriggerMode', 'SetValue', 'RW', PySpin.TriggerMode_Off)
    cam_node_cmd(cam, 'TriggerSource', 'SetValue', 'RW', PySpin.TriggerSource_Software)
    cam_node_cmd(cam, 'TriggerSelector', 'SetValue', 'RW', PySpin.TriggerSelector_FrameStart)
    cam_node_cmd(cam, 'TriggerMode', 'SetValue', 'RW', PySpin.TriggerMode_On)


def set_cam_buffer(cam):
    # PySpin.StreamBufferHandlingMode_NewestFirst
    # PySpin.StreamBufferHandlingMode_NewestOnly
    # PySpin.StreamBufferHandlingMode_OldestFirstOverwrite
    cam.TLStream.StreamBufferHandlingMode.SetValue(PySpin.StreamBufferHandlingMode_NewestOnly)


def disable_auto_exp(cam):
    print('Disabling Auto Exposure')
    cam_node_cmd(cam, 'ExposureAuto', 'SetValue', 'RW', PySpin.ExposureAuto_Off)


def disable_auto_gain(cam):
    print('Disabling Auto Gain')
    cam_node_cmd(cam, 'GainAuto', 'SetValue', 'RW', PySpin.GainAuto_Off)


def configure_trigger(cam, trigger):
    """
    This function configures the camera to use a trigger. First, trigger mode is
    set to off in order to select the trigger source. Once the trigger source
    has been selected, trigger mode is then enabled, which has the camera
    capture only a single image upon the execution of the chosen trigger.
     :param cam: Camera to configure trigger for.
     :type cam: CameraPtr
     :return: True if successful, False otherwise.
     :rtype: bool
    """
    result = True

    print('*** CONFIGURING TRIGGER ***\n')

    print(
        'Note that if the application / user software triggers faster than frame time, the trigger may be dropped / skipped by the camera.\n')
    print(
        'If several frames are needed per trigger, a more reliable alternative for such case, is to use the multi-frame mode.\n\n')

    if trigger == "SOFTWARE":
        print('Software trigger chosen ...')
    elif trigger == "HARDWARE":
        print('Hardware trigger chose ...')

    try:
        # Ensure trigger mode off
        # The trigger must be disabled in order to configure whether the source
        # is software or hardware.
        nodemap = cam.GetNodeMap()
        node_trigger_mode = PySpin.CEnumerationPtr(nodemap.GetNode('TriggerMode'))
        if not PySpin.IsAvailable(node_trigger_mode) or not PySpin.IsReadable(node_trigger_mode):
            print('Unable to disable trigger mode (node retrieval). Aborting...')
            return False

        node_trigger_mode_off = node_trigger_mode.GetEntryByName('Off')
        if not PySpin.IsAvailable(node_trigger_mode_off) or not PySpin.IsReadable(node_trigger_mode_off):
            print('Unable to disable trigger mode (enum entry retrieval). Aborting...')
            return False

        node_trigger_mode.SetIntValue(node_trigger_mode_off.GetValue())

        print('Trigger mode disabled...')

        # Set TriggerSelector to FrameStart
        # For this example, the trigger selector should be set to frame start.
        # This is the default for most cameras.
        node_trigger_selector = PySpin.CEnumerationPtr(nodemap.GetNode('TriggerSelector'))
        if not PySpin.IsAvailable(node_trigger_selector) or not PySpin.IsWritable(node_trigger_selector):
            print('Unable to get trigger selector (node retrieval). Aborting...')
            return False

        node_trigger_selector_framestart = node_trigger_selector.GetEntryByName('FrameStart')
        if not PySpin.IsAvailable(node_trigger_selector_framestart) or not PySpin.IsReadable(
                node_trigger_selector_framestart):
            print('Unable to set trigger selector (enum entry retrieval). Aborting...')
            return False
        node_trigger_selector.SetIntValue(node_trigger_selector_framestart.GetValue())

        print('Trigger selector set to frame start...')

        # Select trigger source
        # The trigger source must be set to hardware or software while trigger
        # mode is off.
        node_trigger_source = PySpin.CEnumerationPtr(nodemap.GetNode('TriggerSource'))
        if not PySpin.IsAvailable(node_trigger_source) or not PySpin.IsWritable(node_trigger_source):
            print('Unable to get trigger source (node retrieval). Aborting...')
            return False

        if trigger == "SOFTWARE":
            node_trigger_source_software = node_trigger_source.GetEntryByName('Software')
            if not PySpin.IsAvailable(node_trigger_source_software) or not PySpin.IsReadable(
                    node_trigger_source_software):
                print('Unable to set trigger source (enum entry retrieval). Aborting...')
                return False
            node_trigger_source.SetIntValue(node_trigger_source_software.GetValue())
            print('Trigger source set to software...')

        elif trigger == "HARDWARE":
            node_trigger_source_hardware = node_trigger_source.GetEntryByName('Line0')
            if not PySpin.IsAvailable(node_trigger_source_hardware) or not PySpin.IsReadable(
                    node_trigger_source_hardware):
                print('Unable to set trigger source (enum entry retrieval). Aborting...')
                return False
            node_trigger_source.SetIntValue(node_trigger_source_hardware.GetValue())
            print('Trigger source set to hardware...')

        # Turn trigger mode on
        # Once the appropriate trigger source has been set, turn trigger mode
        # on in order to retrieve images using the trigger.
        node_trigger_mode_on = node_trigger_mode.GetEntryByName('On')
        if not PySpin.IsAvailable(node_trigger_mode_on) or not PySpin.IsReadable(node_trigger_mode_on):
            print('Unable to enable trigger mode (enum entry retrieval). Aborting...')
            return False

        node_trigger_mode.SetIntValue(node_trigger_mode_on.GetValue())
        print('Trigger mode turned back on...')

    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        return False

    return result


def grab_next_image_by_trigger(nodemap, cam, trigger):
    """
    This function acquires an image by executing the trigger node.
    :param cam: Camera to acquire images from.
    :param nodemap: Device nodemap.
    :type cam: CameraPtr
    :type nodemap: INodeMap
    :return: True if successful, False otherwise.
    :rtype: bool
    """
    try:
        result = True
        # Use trigger to capture image
        # The software trigger only feigns being executed by the Enter key;
        # what might not be immediately apparent is that there is not a
        # continuous stream of images being captured; in other examples that
        # acquire images, the camera captures a continuous stream of images.
        # When an image is retrieved, it is plucked from the stream.

        if trigger == "SOFTWARE":
            # Get user input
            # input('Press the Enter key to initiate software trigger.')

            # Execute software trigger
            node_softwaretrigger_cmd = PySpin.CCommandPtr(nodemap.GetNode('TriggerSoftware'))
            if not PySpin.IsAvailable(node_softwaretrigger_cmd) or not PySpin.IsWritable(node_softwaretrigger_cmd):
                print('Unable to execute trigger. Aborting...')
                return False

            node_softwaretrigger_cmd.Execute()


        elif trigger == "HARDWARE":
            print('Use the hardware to trigger image acquisition.')

    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        return False

    return result
