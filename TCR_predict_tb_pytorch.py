#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 15:16:54 2021

@author: zrollins
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch


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

#Residues
#CAVNGDDYKLSF
#CASRKRDDSEQYF

resi1 = ['A','V','N','G','D','D','Y','K','L','S','F',
        'A','S','R','K','R','D','D','S','E','Q','Y','F',
        'Y','L','S','G','A','N','L','N','L', 3249]

#CAVSDNARLMF
#CASSPFGGGNEQFF

resi2 = ['A','V','S','D','N','A','R','L','M','F','F',
        'A','S','S','P','F','G','G','G','N','E','Q','F',
        'Y','L','S','G','A','N','L','N','L', 1988]

#CAYRISAYDKVIF
#CASSQTGGADTDTQYF

resi3 = ['A','Y','R','I','S','A','Y','D','K','V','I',
        'A','S','S','Q','T','G','G','A','D','T','D','T',
        'Y','L','S','G','A','N','L','N','L', 8639]


data_labels = {'TCR1':resi1, 'TCR2':resi2, 'TCR3':resi3}
dataset = pd.DataFrame(data_labels)
dataset= dataset.T

dataset = dataset.replace(['A','R','N','D','C','Q','E','G','H','I','L',
                         'K','M','F','P','S','T','W','Y','V']
                       ,['1','2','3','4','5','6','7','8','9','10',
                       '11','12','13','14','15','16','17','18','19','20'])

print(dataset)
mutants, xy = dataset.shape
resi = 32
bs = xy - resi
neurons=2

train = dataset.iloc[:3,:]
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


train_set = torch.utils.data.DataLoader(train_data, batch_size=2, shuffle=True)
test_set = torch.utils.data.DataLoader(test_data, batch_size=1, shuffle=False)


class NeuralNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1 = nn.Linear(resi, neurons)
        self.fc2 = nn.Linear(neurons, neurons)
        self.fc3 = nn.Linear(neurons, neurons)
        self.fc4 = nn.Linear(neurons, neurons)
        self.fc5 = nn.Linear(neurons, bs)
        
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        x = F.relu(self.fc4(x))
        x = self.fc5(x)
        return x

net = NeuralNet()
optimizer = optim.Adam(net.parameters(), lr=0.001)
train_losses = []
train_counter = []

for epoch in range(50): # 10 full passes over the data
    for data in train_set:  # `data` is a batch of data
        X, y = data  # X is the batch of features, y is the batch of targets.
        net.zero_grad()  # sets gradients to 0 before loss calc. You will do this likely every step.
        output = net(X)# pass in the batch inputs
        #output = output[:, bs]
        loss = F.mse_loss(output, y)  # calc and grab the loss value
        loss.backward()  # apply this loss backwards thru the network's parameters
        optimizer.step()  # attempt to optimize weights to account for loss/gradients
    train_losses.append(loss.item())
    train_counter.append(epoch)
    print(loss)  # print loss. We hope loss (a measure of wrong-ness) declines! 
    
    
test_loss = 0    
test_losses = []
prediction = np.empty([0, bs])

with torch.no_grad(): #Calculate accuracy of each test mutant
    for data in test_set:
        X, y = data
        output = net(X)
        output = torch.squeeze(output,1)
        test_loss = F.mse_loss(output, y)
        test_losses.append(test_loss.item())
        i=0
        correct = 0
        total = 0
        output = output.numpy()
        print(output)
        prediction = np.append(prediction, output)
        y = y.numpy()
        for pred, act in list(zip(output, y)):
            #print(output, y)
            if abs((pred - act)/act) <= 0.05:
                correct += 1
            total += 1
        print("Accuracy: ", round(correct/total, 3))
      


#prd = pd.DataFrame(prediction, index=['TCR1','TCR2','TCR3'])
#prd.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/fnc_prd.csv', header=None)
 
fig = plt.figure()
plt.plot(train_counter, train_losses, color='blue')
plt.scatter(epoch, test_losses, color='red')
plt.legend(['Train Loss', 'Test Loss'])
plt.ylabel('Batch MSE Loss')
plt.xlabel('Epochs')

