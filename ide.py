
	"""
	#rounding off
	for i in range(3):
		if str(FirstEigenvalues[i])[1] == '.':
			print round(FirstEigenvalues[i],3)
			if int(str(round(FirstEigenvalues[i],3))[-1]) >= 5: #not human readable code
				print "yay"
			else:
				print 'meh1'
		if str(FirstEigenvalues[i])[2] == '.':
			print round(FirstEigenvalues[i],2)
			if int(str(round(FirstEigenvalues[i],3))[-1]) >= 5:
				print "yay2"
			else:
				print 'meh2'
		print '.----'


 
		val = round(FirstEigenvalues[i],3)
		#print val
		#print len(str(val))#if str(val)[4] >= 5:
		#	print "yay"
	#		print str(val)[4]




	"""