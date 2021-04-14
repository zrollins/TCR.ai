#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 13:30:17 2021

@author: zrollins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from math import isclose
import torchvision
from torchvision import transforms, datasets
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim



dataset = pd.read_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/hmapx_pmhc_master_bio.csv', index_col = 0, header=None)
mutants, xy = dataset.shape
resi = 140
bonds = xy - resi
neurons=128

train = dataset.iloc[3:,:]
test = pd.concat([dataset, train], axis=0)
test = test.drop_duplicates(keep=False)

#FUTURE REFERNCE: Consider Embedding residues as n-d arrays

print(train.shape)
print(test.shape)

class Feeder():
    
    def __init__(self,data):
        self.data = torch.FloatTensor(data.values.astype('float'))
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self,index):
        target = self.data[index][resi:]
        data_val =self.data[index] [:resi]
        return data_val, target
    
train_data = Feeder(train)
test_data = Feeder(test)



device = "cuda" if torch.cuda.is_available() else "cpu"
kwargs = {'num_workers': 1, 'pin_memory': True} if device=='cuda' else {}


train_set = torch.utils.data.DataLoader(train_data, batch_size=7, shuffle=True)
test_set = torch.utils.data.DataLoader(test_data, batch_size=3, shuffle=False)


class NeuralNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1 = nn.Linear(resi, neurons)
        self.fc2 = nn.Linear(neurons, neurons)
        self.fc8 = nn.Linear(neurons, bonds)
        self.fc9 = nn.ReLU()
        
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc8(x)
        x = self.fc9(x)
        return x

net = NeuralNet()
optimizer = optim.Adam(net.parameters(), lr=0.001)
train_losses = []
train_counter = []

for epoch in range(25): # 10 full passes over the data
    for data in train_set:  # `data` is a batch of data
        X, y = data  # X is the batch of features, y is the batch of targets.
        net.zero_grad()  # sets gradients to 0 before loss calc. You will do this likely every step.
        output = net(X)# pass in the batch inputs
        #output = output[:, bonds]
        loss = F.mse_loss(output, y)  # calc and grab the loss value
        loss.backward()  # apply this loss backwards thru the network's parameters
        optimizer.step()  # attempt to optimize weights to account for loss/gradients
    train_losses.append(loss.item())
    train_counter.append(epoch)
    print(loss)  # print loss. We hope loss (a measure of wrong-ness) declines! 
    
test_loss = 0    
test_losses = []
prediction = np.empty([0, bonds])

with torch.no_grad(): #Calculate accuracy of each test mutant
    for data in test_set:
        X, y = data
        output = net(X)
        output = torch.squeeze(output,1)
        test_loss = F.mse_loss(output, y)
        print(output)
        print(y)
        test_losses.append(test_loss.item())
        i=0
        correct = 0
        total = 0
        output = output.numpy()
        print(output)
        prediction = np.append(prediction, output, axis=0)
        y = y.numpy()
        print(output.shape)
        print(y.shape)
        print(prediction.shape)
      


prd = pd.DataFrame(prediction, index=['L1','MART1','GVA'])
prd.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/hmapx_pmhc_mstr_prd_2x128_bio.csv', header=None)
 
fig = plt.figure()
plt.plot(train_counter, train_losses, color='blue')
plt.scatter(epoch, test_losses, color='red')
plt.legend(['Train Loss', 'Test Loss'])
plt.ylabel('Batch MSE Loss')
plt.xlabel('Epochs')