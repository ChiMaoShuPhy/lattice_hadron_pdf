from lxml import etree
with open('try_in.xml', 'r') as xml_file:
    xml_string=xml_file.read().replace('\n', '')
root = etree.fromstring(xml_string)
tree = etree.ElementTree(root)
for e in root.iter():
    if (e.tag=='SinkType'):
        print tree.getpath(e)