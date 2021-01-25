from pycaret.datasets import get_data
juice = get_data('juice')
print(juice.Purchase.unique())

#exp_name = setup(data = juice,  target = 'Purchase')
'''lr = create_model('lr')
plot_model(lr, plot = 'auc')'''