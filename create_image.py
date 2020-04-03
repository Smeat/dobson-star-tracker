#!/usr/bin/env python3

from argparse import ArgumentParser, RawTextHelpFormatter
from enum import Enum
from typing import List
from collections import defaultdict

import random
import operator
import math
import csv
import re
import struct

INDEX_ENTRY_SIZE = 12

MAGIC_ID = 0x11223344
VERSION = 0


def to_sector(pos):
	return math.ceil(pos/512)

def to_pos(sector):
	return sector * 512

def next_sector_pos(pos):
	return pos + (512 - (pos % 512))

class SkyObjectType(Enum):
	UNDEFINED = 0
	STAR = 1
	GALAXY = 2
	CLUSTER_GLOBULAR = 3
	CLUSTER_OPEN = 4
	NEBULA = 5

class ImageHeader:
	SIZE = 12
	magic_id: int
	version: int
	type_begin: int

	def __init__(self, version, type_begin):
		self.magic_id = MAGIC_ID
		self.version = version
		self.type_begin = type_begin

	def pack(self):
		return struct.pack('>III', self.magic_id, self.version, self.type_begin)

class ImageTypeHeader:
	SIZE = 9
	sky_type: SkyObjectType
	data_begin: int
	name_begin: int

	def __init__(self, sky_type, data_begin, name_begin):
		self.sky_type = sky_type
		self.data_begin = data_begin
		self.name_begin = name_begin

	def pack(self):
		return struct.pack('>BII', self.sky_type, self.data_begin, self.name_begin)

class ImageEntry:
	SIZE = 16
	ID: int
	RA: int
	dec: int
	mag: int
	flags: int

	def __init__(self, RA, dec, mag, flags, ID = None):
		if ID is None:
			ID = random.randint(0, 2**31)
		self.ID = ID
		self.RA = RA
		self.dec = dec
		self.mag = mag
		self.flags = flags

	def __str__(self):
		return "{}".format(self.mag)

	def pack(self):
		return struct.pack('>IiihH', self.ID, self.RA, self.dec, self.mag, self.flags)

class ImageNameEntry:
	SIZE = 40
	shortname: str
	longname: str

	def __init__(self, shortname, longname):
		self.shortname = shortname
		self.longname = longname

	def pack(self):
		return struct.pack('8s32s', self.shortname, self.longname)


class SkyObject:
	objectType: SkyObjectType
	RA: int
	DEC: int
	identifier: str
	name: str
	mag: int

	def __init__(self, RA, DEC, identifier, name, mag, objectType = 0 ):
		self.objectType = objectType
		self.RA = RA
		self.DEC = DEC
		self.identifier = identifier
		self.name = name
		self.mag = mag


def kstars_binary_read_header(catalog):
	header_text = catalog.read(124)
	print("Header is: {}".format(header_text))
	endianess = catalog.read(2)
	if(endianess != b'SK'):
		print("Byteswapping not implemented!")
		return
	version = struct.unpack('b', catalog.read(1))[0]
	print("Version: {}".format(version))
	num_fields = struct.unpack('h', catalog.read(2))[0]
	print("Number of fields: {}".format(num_fields))
	for i in range(0, num_fields):
		field = struct.unpack('10sbBi', catalog.read(16))
		print("Field: {}".format(field))
	trixels = struct.unpack('i', catalog.read(4))[0]
	print("Number of trixels: {}".format(trixels))
	print("header done! file pos: {}".format(catalog.tell()))
	return (catalog.tell(), trixels)


def parse_kstars_binary(filename, starname_file = "/usr/share/kstars/starnames.dat") -> List[SkyObject]:
	sky_objects = []
	with open(filename, "rb") as catalog, open(starname_file, "rb") as starnames:
		header_size, trixels = kstars_binary_read_header(catalog)
		header_size_name = kstars_binary_read_header(starnames)[0]
		# header done
		starnames.seek(header_size_name + 1 * INDEX_ENTRY_SIZE)
		for t in range(0, trixels):
			catalog.seek(header_size + t * INDEX_ENTRY_SIZE) # seek to a specific trixel
			f_trix = struct.unpack('i', catalog.read(4))[0]
			print("Read trixel: {}".format(f_trix))
			offset = struct.unpack('i', catalog.read(4))[0]
			nrec = struct.unpack('i', catalog.read(4))[0]
			print("offset: {}, nrec: {}".format(offset, nrec))
			catalog.seek(offset)
			for i in range(0, nrec):
				offset = catalog.tell()
				star_data = struct.unpack('iiiiiihhhbb', catalog.read(32))
				print("Star data: {}".format(star_data))
				if star_data[9] & 0x01:
					print("has name!")
					shortname = struct.unpack('8s', starnames.read(8))[0]
					longname = struct.unpack('32s', starnames.read(32))[0]
					print("Shortname: {} Longname: {}".format(shortname, longname))
				sky_objects.append(SkyObject(RA = star_data[0], DEC = star_data[1], mag = star_data[6], identifier = shortname, name = longname))
	return sky_objects


def parse_kstars_catalog(filename):
	with open(filename) as catalog:
		# find field name
		pattern = re.compile("^#.* RA .*DEC.*")
		for row in catalog:
			match = pattern.match(row)
			if match is not None:
				fields = match.group().split(" ")
				fields = list(filter(None, fields))[1:]
				break

		csv_reader = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter=' ', fieldnames=fields)
		for row in csv_reader:
			print("Row: {}".format(row))

def write_image(outputfilename, objects: List[SkyObject]):
	objects.sort(key=operator.attrgetter('mag'))
	with open(outputfilename, "wb") as outfile:
		image_header = ImageHeader(0, to_sector(ImageHeader.SIZE))
		outfile.write(image_header.pack())
		# gather types
		types = defaultdict(lambda: [])
		names = []
		for o in objects:
			entry = ImageEntry(RA=o.RA, dec = o.DEC, mag=o.mag, flags= 1)
			name_entry = ImageNameEntry(o.identifier, o.name)
			types[o.objectType].append([entry, name_entry])
			names.append(name_entry)
		# they should already be sorted by mag

		# total size for all type entries
		type_entry_size = len(types.keys()) * ImageTypeHeader.SIZE
		type_entry_sectors = math.ceil(type_entry_size / 512)
		data_entry_size = len(objects) * ImageEntry.SIZE
		data_entry_sectors = math.ceil(data_entry_size / 512)

		type_entry_begin = to_pos(image_header.type_begin)
		data_entry_begin = next_sector_pos(type_entry_begin + type_entry_size)
		name_entry_begin = next_sector_pos(data_entry_begin + data_entry_size)

		print("Type entry begins at {} data {} name {}".format(type_entry_begin, data_entry_begin, name_entry_begin))

		# write type list
		outfile.seek(type_entry_begin)
		next_data_entry_begin = data_entry_begin
		next_name_entry_begin = name_entry_begin
		for t, e in types.items():
			type_list = ImageTypeHeader(t, next_data_entry_begin >> 9, next_name_entry_begin >> 9)
			print("Writing type header at pos: {}".format(outfile.tell()))
			outfile.write(type_list.pack())
			outfile.seek(next_data_entry_begin)
			for entry, name in e:
				outfile.write(entry.pack())
			outfile.seek(next_sector_pos(outfile.tell())) # seek to begin of next sector
			next_data_entry_begin = outfile.tell()
			outfile.seek(next_name_entry_begin)
			for entry, name in e:
				outfile.write(name.pack())
			outfile.seek(next_sector_pos(outfile.tell())) # seek to begin of next sector
			next_name_entry_begin = outfile.tell()


		print("type entry size {} takes {} sectors. data entry size {} takes {} sectors".format(type_entry_size, type_entry_sectors, data_entry_size, data_entry_sectors))



if __name__ == "__main__":
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("-i", "--input", dest="input", help="Path to input file", type=str, required=True)
	args = parser.parse_args()

	kstars_data = parse_kstars_binary(args.input)
	write_image("/tmp/test", kstars_data)
