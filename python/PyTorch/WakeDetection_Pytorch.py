#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from torchvision import datasets, transforms
from torch.utils.data import Dataset,DataLoader,random_split
batch_size = 16


#root directory
datadir = os.path.expanduser("/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4_figs_thr70")

# #transformations, options
# train_transforms = transforms.Compose([transforms.Resize((224,224)),
#                                        transforms.ToTensor(),                                
#                                        torchvision.transforms.Normalize(
#                                            mean=[0.485, 0.456, 0.406],
#                                            std=[0.229, 0.224, 0.225],
#     ),
#                                        ])
# test_transforms = transforms.Compose([transforms.Resize((224,224)),
#                                       transforms.ToTensor(),
#                                       torchvision.transforms.Normalize(
#                                           mean=[0.485, 0.456, 0.406],
#                                           std=[0.229, 0.224, 0.225],
#     ),
#                                       ])

#transformations, simpler
transforms = transforms.Compose([transforms.Resize((224,224)),
                                        transforms.ToTensor()])


#datasets
dataset = datasets.ImageFolder(datadir, transform=transforms)
# dataset = datasets.ImageFolder(datadir)

# Define the proportion or number of items in each set
train_size = int(0.8 * len(dataset))
test_size = len(dataset) - train_size

# Randomly split the dataset into train and test datasets
training_data, test_data = random_split(dataset, [train_size, test_size])


# Create data loaders.
train_dataloader = DataLoader(training_data, batch_size=batch_size, shuffle=True)
test_dataloader = DataLoader(test_data, batch_size=batch_size, shuffle=False)

for X, y in train_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break
    
for X, y in test_dataloader:
    print(f"Shape of X [N, C, H, W]: {X.shape}")
    print(f"Shape of y: {y.shape} {y.dtype}")
    break    


# In[2]:


import torch
from torch import nn

# Get cpu or gpu device for training.
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")

# Define model
class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(224*224*3, 32),
            nn.ReLU(),
#             nn.Linear(512, 512),
#             nn.ReLU(),
            nn.Linear(32, 2)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits

model = NeuralNetwork().to(device)
print(model)


# In[3]:


loss_fn = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=1e-3)


# In[4]:


def train(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    model.train()
    for batch, (X, y) in enumerate(dataloader):
        X, y = X.to(device), y.to(device)

        # Compute prediction error
        pred = model(X)
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch % 100 == 0:
            loss, current = loss.item(), batch * len(X)
            print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")


# In[5]:


def test(dataloader, model, loss_fn):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    model.eval()
    test_loss, correct = 0, 0
    with torch.no_grad():
        for X, y in dataloader:
            X, y = X.to(device), y.to(device)
            pred = model(X)
            test_loss += loss_fn(pred, y).item()
            correct += (pred.argmax(1) == y).type(torch.float).sum().item()
    test_loss /= num_batches
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")


# In[6]:


epochs = 5
for t in range(epochs):
    print(f"Epoch {t+1}\n-------------------------------")
    train(train_dataloader, model, loss_fn, optimizer)
    test(test_dataloader, model, loss_fn)
print("Done!")


# In[ ]:




