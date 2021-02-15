import numpy as np
import torch
from scipy.signal import resample
from PIL import Image, ImageDraw
from frospy.util.Convolutional.ConvNetwork import Net
from torchvision import transforms
import torch.nn.functional as F
from frospy import data as frospydata


class ConvNet:
    def __init__(self):
        # load torch model
        model_path = "%s/cnn/model_dict.pth" % frospydata.__path__[0]
        model = Net()
        model.load_state_dict(torch.load(model_path))
        model.eval()
        self.model = model

    def Predict(self, data):
        line = []

        # resample data to size used for training and normalise
        d = resample(data, 15)
        d = (d - np.min(d)) / (np.max(d) - np.min(d))*15

        # create an empty black square image
        img = Image.new('RGB', (len(d)-1, len(d)+1), 'black')
        # create line coordinates from data
        for i in range(0, len(d)-1):
            line.append(i)
            line.append(d[i])

        # plot (or draw) line onto black square
        draw = ImageDraw.Draw(img)
        draw.line(line, fill='white', width=1)
        img = img.transpose(Image.ROTATE_180)

        trans = transforms.ToTensor()
        img_tensor = trans(img).unsqueeze(0)

        # predict class of pick based on highest class probability
        output = self.model(img_tensor)
        index = output.data.cpu().numpy().argmax()
        prob = F.softmax(output)
        prob = prob.detach().numpy()
        prob = np.unique(prob)
        return index, prob
