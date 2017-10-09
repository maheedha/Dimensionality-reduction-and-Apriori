
# coding: utf-8

# In[2]:

import numpy as np
import itertools
import collections
from collections import defaultdict
import re

# Apriori Algorithm Implementation

results = []
minsup_percent = 50
len_list = []
dict_list = []
with open('associationruletestdata.txt') as inputfile:
    for line in inputfile:
        each_line = line.strip().split('\t')
        line_data = []
        count = 1
        for value in each_line:
            if(count < 10):
                line_data.append('G' + str(0) + str(count) + '_' + value)
            else:
                line_data.append('G' + str(count) + '_' + value)
            count += 1
        results.append(line_data)
minsup = int((len(results[0]) - 1) * minsup_percent /100)
dataset = np.asarray(results)
unique_itemsets = np.unique(dataset)
#  k = 1
unique, counts = np.unique(results, return_counts=True)
unique_dict = dict(zip(unique, counts))
len_list.append(1)
dict_list.append(unique_dict)
print('Support is set to be ' +str(minsup_percent) + '%')
support = (counts>=minsup).sum()

print('number of length-1 frequent itemsets:' +str(support))
unique_dict = {k: v for k,v in unique_dict.items() if v >= minsup}
unique_keys = unique_dict.keys()
# k = 2
two_sets = []
two_sets = itertools.combinations(unique_keys,2)
two_sets = [list(t) for t in two_sets]
two_sets = np.array(two_sets)
unique_dict_k = {}
for patient in results:
    for comb in two_sets:
        if(set(comb).issubset(set(patient))):
            temp = ','.join(comb)
            if temp in unique_dict_k:
                unique_dict_k[temp] += 1
            else:
                unique_dict_k[temp] = 1
support = 0
for k,v in unique_dict_k.items():
    if v >= minsup:
        support += 1
print('number of length-2 frequent itemsets:' +str(support))
unique_dict_k = {k: v for k,v in unique_dict_k.items() if v >= minsup}    
len_list.append(2)
dict_list.append(unique_dict_k)

iteration = 2
# from k = 3 to end
index = 1
while(support > 0):
    iteration += 1
    support = 0
    unique_keys = []
    temp_list = []
    for k,v in unique_dict_k.items():
        unique_keys.append(list(n for n in k.split(',')))
    for i in range(0, len(unique_keys)-1):
        for j in range(i+1, len(unique_keys)):
            key1 = unique_keys[i]
            key2 = unique_keys[j]
            flag = True
            for z in range(0,index):
                if key1[z] != key2[z]:
                    flag = False
            if(flag):
                t = list(set(key1) | set(key2))
                t.sort()
                temp_list.append(t)
    unique_dict_k = {}
    for patient in results:
        for comb in temp_list:
            if(set(comb).issubset(set(patient))):
                temp = ','.join(comb)
                if temp in unique_dict_k:
                    unique_dict_k[temp] += 1
                else:
                    unique_dict_k[temp] = 1
    index += 1
    for k,v in unique_dict_k.items():
        if v >= minsup:
            support += 1
    print('number of length-' + str(iteration) +' frequent itemsets:' +str(support))
    unique_dict_k = {k: v for k,v in unique_dict_k.items() if v >= minsup}
    unique_dict_k = collections.OrderedDict(sorted(unique_dict_k.items()))
#     print(unique_dict_k)
    if(support > 0):
        len_list.append(iteration)
        dict_list.append(unique_dict_k)
        
#  Rule generation
        
confup_percent = 70
confidence = int((len(results[0]) - 1) * confup_percent /100)
intr_list = []
conf_list = []
unique_keys = []
total_iterations = len(len_list)-1
# print(len(len_list))
while(total_iterations > 0):
    curr_dict = dict_list[total_iterations]
    for k,v in curr_dict.items():
            temp = []
            temp.append([n for n in k.split(',')])
            for i in range(0, len(temp[0])):
                temp_dict = defaultdict(list)
                for j in range(0, len(temp[0])):
                    if(temp[0][i] != temp[0][j]):
                        temp_dict[temp[0][i]].append(temp[0][j])
                val = temp_dict[temp[0][i]]
                val.sort()
                flag1 = False
                val = ','.join(val)
                k_set = dict_list[total_iterations-1]
                if val in k_set and ((v / k_set[val]) * 100 >= confidence):
                    intr_list.append({val:temp[0][i]})
                    flag1 = True
            if flag1:
                for i in range(2, total_iterations+1):
                    l1 = k.split(',')
                    st = ','.join(l1)
                    perm = itertools.combinations(l1,i)
                    k_set = dict_list[len(len_list)-i-1]
                    perm = [list(t) for t in perm]
                    for i in range(0, len(perm)):
                        l2 = [x for x in l1 if x not in perm[i]]
                        str1 = ','.join(perm[i])
                        str2 = ','.join(l2)
                        if str2 in k_set:
                            if(float((v / k_set[str2]) * 100 >= confidence)):
                                intr_list.append({str2:str1})
    total_iterations -= 1

# template functions
def any_fn(rules):
    res = []
    if rules[0].strip() in "RULE":
        for dict1 in intr_list:
            for key, value in dict1.items():
                for i in range(2, len(rules)):
                    if rules[i].strip() in key or  rules[i].strip() in value:
                        res.append({key:value})
    elif rules[0].strip() in "BODY":
        for dict1 in intr_list:
            for key,value in dict1.items():
                for i in range(2, len(rules)):
                    if rules[i].strip() in key:
                        res.append({key:value})
    elif rules[0].strip() in "HEAD":
        for dict1 in intr_list:
            for key,value in dict1.items():
                for i in range(2, len(rules)):
                    if rules[i].strip() in value:
                        res.append({key:value})
    return res

def none_fn(rules):
    res = []
    if rules[0].strip() in "RULE":
            for dict1 in intr_list:
                for key, value in dict1.items():
                    flag = False
                    for i in range(2, len(rules)):
                        if rules[i].strip() in key or rules[i].strip() in value:
                            flag = True
                    if(flag):
                        continue
                    res.append({key:value})
    elif rules[0].strip() in "BODY":
        for dict1 in intr_list:
            for key,value in dict1.items():
                flag = False
                for i in range(2, len(rules)):
                    if rules[i].strip() in key:
                        flag = True
                if(flag):
                    continue
                res.append({key:value})
    elif rules[0].strip() in "HEAD":
        for dict1 in intr_list:
            for key,value in dict1.items():
                flag = False
                for i in range(2, len(rules)):
                    if rules[i].strip() in value:
                        flag = True
                if(flag):
                    continue
                res.append({key:value})
    return res

def other_fn(rules):
    res = []
    if len(rules) == 3:
        res = any_fn(rules)
    elif rules[0].strip() in "RULE":
        for dict1 in intr_list:
            count = 0
            for key, value in dict1.items():
                for i in range(2, len(rules)):
                    if bool(rules[i].strip() in key) != bool(rules[i].strip() in value):
                        count += 1
                if count == 1:
                    res.append(dict1)
    elif rules[0].strip() in "BODY":
        for dict1 in intr_list:
            count = 0
            for key,value in dict1.items():
                for i in range(2, len(rules)):
                    if bool(rules[i].strip() in key):
                        count += 1
                if count == 1:
                    res.append(dict1)
    elif rules[0].strip() in "HEAD":
        for dict1 in intr_list:
            count = 0
            for key,value in dict1.items():
                for i in range(2, len(rules)):
                    if bool(rules[i].strip() in value):
                        count += 1
                if count == 1:
                    res.append(dict1)
    return res

def template1_caller(rules):
    res = []
    if rules[1] == ' ANY':
        res = any_fn(rules)
    elif rules[1] == ' NONE':
        res = none_fn(rules)
    elif rules[1]:
        res = other_fn(rules)
    return res

# template2
def template2_caller(rules,num):
    res = []
    if rules == ' RULE':
        for dict_l in intr_list:
                for k,v in dict_l.items():
                    keys = k.split(',')
                    values = v.split(',')
                    dict_len = len(keys) + len(values)
                    if dict_len >= int(num) :
                        res.append(dict_l)
    elif rules == ' BODY':
        for dict_l in intr_list:
                for k,v in dict_l.items():
                    keys = k.split(',')
                    values = v.split(',')
                    dict_len = len(keys)
                    if dict_len >= int(num) :
                        res.append(dict_l)
    elif rules == ' HEAD':
        for dict_l in intr_list:
                for k,v in dict_l.items():
                    keys = k.split(',')
                    values = v.split(',')
                    dict_len = len(values)
                    if dict_len >= int(num) :
                        res.append(dict_l)
    return res

def dic_union(res1,res2) :
    res_union = []
    for dict1 in res1:
        for k,v in dict1.items():
            if(dict1 not in res_union):
                res_union.append({k:v})
    for dict2 in res2:
        for k,v in dict2.items():
            if(dict2 not in res_union):
                res_union.append({k:v})
    return res_union

def dic_intersec(res1,res2) :
    res_intersection = []
    for dict1 in res1:
        for dict2 in res2:
            if dict1 == dict2:
                res_intersection.append(dict1)
    return(res_intersection)

# Querying
print()
# template1
with open('template1.txt') as inputfile:
    for line in inputfile:
        line1 = line.replace("asso_rule.template1", "")
        resul, val = line1.split("=")
        val = re.sub("[^a-z0-9_A-Z, ]+", "", val)
        rules = val.split(',')
        res = template1_caller(rules) 
        print(line +" : " + str(res) + "and length is " +str(len(res)))
        print()  

# template2
with open('template2.txt') as inputfile:
    for line in inputfile:
        line1 = line.replace("asso_rule.template2", "")
        res, val = line1.split("=")
        val = re.sub("[^a-z0-9_A-Z, ]+", "", val)
        rules = val.split(',')
        left=rules[0]
        right=rules[1]
        res = template2_caller(left,right)
        print(line +" : " + str(res) + "and length is " +str(len(res)))
        print()
        
# template3
with open('template3.txt') as inputfile:
    for line in inputfile:
        line1 = line.replace("asso_rule.template3", "")
        resul, val = line1.split("=")
        val = re.sub("[^a-z0-9_A-Z, ]+", "", val)
        rules = val.split(',')
        res_intersection = []
        res_union = []
        if rules[0].strip() == '1or1':
            a = rules[1:4]
            b = rules[4:]
            res1 = template1_caller(a)
            res2 = template1_caller(b)
            for dict1 in res1:
                for k,v in dict1.items():
                    if(dict1 not in res_union):
                        res_union.append({k:v})
            for dict2 in res2:
                for k,v in dict2.items():
                    if(dict2 not in res_union):
                        res_union.append({k:v})
            print(line +" : " + str(res_union) + "and length is " +str(len(res_union)))
            print()
        elif rules[0].strip() == '1and1':
            a = rules[1:4]
            b = rules[4:]
            res1 = template1_caller(a)
            res2 = template1_caller(b)
            for dict1 in res1:
                for dict2 in res2:
                    if dict1 == dict2:
                        res_intersection.append(dict1)
            print(line +" : " + str(res_intersection) + "and length is " +str(len(res_intersection)))
            print()
        elif rules[0].strip() =='1or2' :
            a = rules[1:4]
            res1 = template1_caller(a)
            res2 = template2_caller(rules[4],rules[5])
            for dict1 in res1:
                for k,v in dict1.items():
                    if(dict1 not in res_union):
                        res_union.append({k:v})
            for dict2 in res2:
                for k,v in dict2.items():
                    if(dict2 not in res_union):
                        res_union.append({k:v})
            print(line +" : " + str(res_union) + "and length is " +str(len(res_union)))
            print()
        elif rules[0].strip() =='1and2' :
            a = rules[1:4]
            res1 = template1_caller(a)
            res2 = template2_caller(rules[4],rules[5])
            res_intersection = dic_intersec(res1,res2)
            print(line +" : " + str(res_intersection) + "and length is " +str(len(res_intersection)))
            print()
        elif rules[0].strip() == '2or2' :
            res1 = template2_caller(rules[1],rules[2])
            res2 = template2_caller(rules[3],rules[4])
            res_union = dic_union(res1,res2)
            print(line +" : " + str(res_union) + "and length is " +str(len(res_union)))
            print()
        elif rules[0].strip() == '2and2' :
            res1 = template2_caller(rules[1],rules[2])
            res2 = template2_caller(rules[3],rules[4])
            res_intersection = dic_intersec(res1,res2)
            print(line +" : " + str(res_intersection) + "and length is " +str(len(res_intersection)))


# In[ ]:



