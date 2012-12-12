
#Store the base pairs
for i in len(str)
	if i = '(':
		push(i)
	if i = ')':
		j = pop()
		addpair(j,i)

#calculate the break number

b=0
for each pair w1,w2 in WTset
	for each pair c1, c2 in Constraintset
		if w1==c1 && w2==c2: continue
		elif w1==c1 && w2!=c2: b++
		elif w1!=c1 && w2==c2: b++ 

		
		elif w1>c1 && w1<c2:
			if w2>c2: b++ 

		elif w2>c1 && w2<c2:
			if w1<c1: b++


return b
