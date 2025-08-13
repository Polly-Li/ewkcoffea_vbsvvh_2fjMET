#put csv json stuffs

def json_to_dict(file_name):
    '''load json file as dict'''
    import json
    with open(file_name, 'r') as f:
        data = json.load(f)
    return data

def save_dict_to_json(data_dict, file_name):
    import json
    from common import create_folder
    """save a (nested) json as json with given file path+name"""
    create_folder(file_name)
    try:
        with open(file_name, 'w') as f:
            json.dump(data_dict, f, indent=4)
        print(f"Dictionary successfully saved to '{file_name}'")
    except Exception as e:
        print(f"Failed to save dictionary to '{file_name}': {e}")
        
def save_array_to_csv(arr, file_name):
    import csv
    """
    Save a 2D array to a CSV file.
    Assumes arr[0] is the header row, arr[1:] are data rows.
    """
    with open(file_name, mode='w', newline='') as f:
        writer = csv.writer(f)
        for row in arr:
            writer.writerow(row)