#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 13:10:23 2024

@author: asus
"""



#%%


# parser
# Run example:
# python WakeDetection_Pytorch_CNN_TranfLearn_term --batch_size=32 --num_workers=0



import argparse

parser = argparse.ArgumentParser(description='efficientnet_b7 wake classification model')
# parser.add_argument('--lr', default=0.1, help='')
parser.add_argument('--batch_size', type=int, default=32, help='')
parser.add_argument('--num_workers', type=int, default=0, help='')
parser.add_argument('--num_epochs', type=int, default=10, help='')

args = parser.parse_args()

# parameters
batch_size =  args.batch_size
print("Batch size = "+ str(batch_size))

num_workers=args.num_workers
print("Num of CPU workers = "+ str(num_workers))

num_epochs=args.num_epochs
print("Num epochs = "+ str(num_epochs))

# batch_size = 32
TRAIN_RATIO = 0.8       #fraction of total dataset that will go to train+validation
VALID_RATIO = 0.9       #fraction of train+validation dataset that will *NOT* go to validation
OUTPUT_DIM = 1          # 2 classes for classification labels
SEED = 1234
pretrained_size = 512

#%%
# Define the data transformations:

    
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


