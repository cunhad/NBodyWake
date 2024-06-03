#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import the modules we need

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.optim.lr_scheduler import _LRScheduler
import torch.utils.data as data
from torch.utils.data import Dataset,DataLoader,random_split

import torchvision.transforms as transforms
import torchvision.datasets as datasets

from sklearn import decomposition
from sklearn import manifold
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
from tqdm.notebook import tqdm, trange
import matplotlib.pyplot as plt
import numpy as np

import copy
import random
import time


# In[2]:


# set the random seeds

SEED = 1234

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
torch.backends.cudnn.deterministic = True


# In[3]:


#data
import os


batch_size = 8
TRAIN_RATIO = 0.8       #fraction of total dataset that will go to train+validation
VALID_RATIO = 0.9       #fraction of train+validation dataset that will *NOT* go to validation
OUTPUT_DIM = 2          # 2 classes for classification labels


#root directory
datadir = os.path.expanduser("/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr70")

pretrained_size = 224
pretrained_means = [0.485, 0.456, 0.406]
pretrained_stds = [0.229, 0.224, 0.225]

train_transforms = transforms.Compose([
                           transforms.Resize(pretrained_size),
                           transforms.RandomRotation(5),
                           transforms.RandomHorizontalFlip(0.5),
                           transforms.RandomCrop(pretrained_size, padding=10),
                           transforms.ToTensor(),
                           transforms.Normalize(mean=pretrained_means,
                                                std=pretrained_stds)
                       ])

test_transforms = transforms.Compose([
                           transforms.Resize(pretrained_size),
                           transforms.ToTensor(),
                           transforms.Normalize(mean=pretrained_means,
                                                std=pretrained_stds)
                       ])



#datasets
# dataset = datasets.ImageFolder(datadir, transform=transforms)
dataset = datasets.ImageFolder(datadir)
print(dataset.class_to_idx)   # this is the label dictionary

# Define the proportion or number of items in each set
train_size = int(TRAIN_RATIO * len(dataset))
test_size = len(dataset) - train_size

# Randomly split the dataset into train and test datasets
train_data, test_data = random_split(dataset, [train_size, test_size])
test_data.dataset.transform = test_transforms
train_data.dataset.transform = train_transforms

# Randomly split the train+validation

n_train_examples = int(len(train_data) * VALID_RATIO)
n_valid_examples = len(train_data) - n_train_examples

train_data, valid_data = data.random_split(train_data,
                                           [n_train_examples, n_valid_examples])

# Create data loaders.
train_dataloader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
test_dataloader = DataLoader(test_data, batch_size=batch_size, shuffle=False)
valid_dataloader = DataLoader(valid_data, batch_size=batch_size, shuffle=False)

for X, y in train_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break
    
for X, y in test_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break  

for X, y in valid_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break  



# In[4]:


# # plot out a few images to ensure the transformations look sensible

# def normalize_image(image):
#     image_min = image.min()
#     image_max = image.max()
#     image.clamp_(min=image_min, max=image_max)
#     image.add_(-image_min).div_(image_max - image_min + 1e-5)
#     return image


# def plot_images(images, labels, classes, normalize=True):

#     n_images = len(images)

#     rows = int(np.sqrt(n_images))
#     cols = int(np.sqrt(n_images))

#     fig = plt.figure(figsize=(10, 10))

#     for i in range(rows*cols):

#         ax = fig.add_subplot(rows, cols, i+1)

#         image = images[i]

#         if normalize:
#             image = normalize_image(image)

#         ax.imshow(image.permute(1, 2, 0).cpu().numpy())
#         ax.set_title(classes[labels[i]])
#         ax.axis('off')


# N_IMAGES = 3

# images, labels = zip(*[(image, label) for image, label in
#                        [train_data[i] for i in range(N_IMAGES)]])

# classes = test_data

# plot_images(images, labels, classes)


# In[5]:


# Defining the Model

# define the VGG base architecture

class VGG(nn.Module):
    def __init__(self, features, output_dim):
        super().__init__()

        self.features = features

        self.avgpool = nn.AdaptiveAvgPool2d(7)

        self.classifier = nn.Sequential(
            nn.Linear(512 * 7 * 7, 4096),
            nn.ReLU(inplace=True),
            nn.Dropout(0.5),
            nn.Linear(4096, 4096),
            nn.ReLU(inplace=True),
            nn.Dropout(0.5),
            nn.Linear(4096, output_dim),
        )

    def forward(self, x):
        x = self.features(x)
        x = self.avgpool(x)
        h = x.view(x.shape[0], -1)
        x = self.classifier(h)
        return x, h

# Below are the configurations for VGG11, VGG13, VGG16 and VGG19. 

vgg11_config = [64, 'M', 128, 'M', 256, 256, 'M', 512, 512, 'M', 512, 512, 'M']

vgg13_config = [64, 64, 'M', 128, 128, 'M', 256, 256, 'M', 512, 512, 'M', 512,
                512, 'M']

vgg16_config = [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 'M', 512, 512, 512,
                'M', 512, 512, 512, 'M']

vgg19_config = [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 256, 'M', 512, 512,
                512, 512, 'M', 512, 512, 512, 512, 'M']

# define a function which takes in a configuration list and returns a nn.Sequential

def get_vgg_layers(config, batch_norm):

    layers = []
    in_channels = 3

    for c in config:
        assert c == 'M' or isinstance(c, int)
        if c == 'M':
            layers += [nn.MaxPool2d(kernel_size=2)]
        else:
            conv2d = nn.Conv2d(in_channels, c, kernel_size=3, padding=1)
            if batch_norm:
                layers += [conv2d, nn.BatchNorm2d(c), nn.ReLU(inplace=True)]
            else:
                layers += [conv2d, nn.ReLU(inplace=True)]
            in_channels = c

    return nn.Sequential(*layers)

# get the features for the VGG11 architecture, with batch normalization.

vgg11_layers = get_vgg_layers(vgg11_config, batch_norm=True)

print(vgg11_layers)

# pass these features to our base VGG module to get our VGG11 model

OUTPUT_DIM = 2

model = VGG(vgg11_layers, OUTPUT_DIM)

print(model)


# In[6]:


# Pre-trained Models

# Let's import a pre-trained VGG11 with batch normalization.

import torchvision.models as models

pretrained_model = models.vgg11_bn(pretrained=True)

print(pretrained_model)

pretrained_model.classifier[-1]


# In[7]:


# define a new final linear layer 

IN_FEATURES = pretrained_model.classifier[-1].in_features

final_fc = nn.Linear(IN_FEATURES, OUTPUT_DIM)

# overwrite the previous linear layer with our new linear layer

pretrained_model.classifier[-1] = final_fc

print(pretrained_model.classifier)


# In[8]:


# load the parameters of the pretrained_model into the model

model.load_state_dict(pretrained_model.state_dict())

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


print(f'The model has {count_parameters(model):,} trainable parameters')


# In[9]:


# # freeze parameters 

# # if we wanted to freeze the features layer then we could do that with

# for parameter in model.features.parameters():
#     parameter.requires_grad = False

# # Freezing all but the last layer in the classifier can be done with:
# for parameter in model.classifier[:-1].parameters():
#     parameter.requires_grad = False



# In[10]:


# Training the Model

START_LR = 1e-7

optimizer = optim.Adam(model.parameters(), lr=START_LR)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

criterion = nn.CrossEntropyLoss()

model = model.to(device)
criterion = criterion.to(device)


# In[]:

    



