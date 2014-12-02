//Bulk code by Thomas
#pragma once

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include <iostream>
#include <vector>

class LayerClassification {

public:

	enum BarrelTypes {
		PXB = 1, TIB = PXB + 3, TOB = TIB + 4, bUNKNOWN = TOB + 6
	};
	enum EndcapTypes {
		PXF = 101, TID = PXF + 2, TEC = TID + 3, eUNKNOWN = TID + 9
	};


	typedef std::vector<LayerClassification> tLayers;

	LayerClassification() {
		//set default
		barrel = false;
		barrelType = BarrelTypes::bUNKNOWN;
		endcapType = EndcapTypes::eUNKNOWN;
		layer = 0;
		side = 0;
		stereo = false;
	}

	LayerClassification(const DetId& detID) {

		//set default
		barrel = false;
		barrelType = BarrelTypes::bUNKNOWN;
		endcapType = EndcapTypes::eUNKNOWN;
		layer = 0;
		side = 0;
		stereo = false;

		//only works for tracker for now
		if (detID.det() != DetId::Detector::Tracker) {
			std::cout << "Error: Bad DetID, only Tracker for now" << std::endl;
			return;
		}

		//copied from Alignment_OfflineValidation_TrackerGeometryCompare.cc
		switch (detID.subdetId()) {

		case PixelSubdetector::PixelBarrel: {
			PXBDetId sid(detID.rawId());
			barrel = true;
			barrelType = BarrelTypes::PXB;
			layer = sid.layer() + barrelType - 1;
			stereo = true; //pixel is 3D
			ladder = sid.ladder();
			break;
		}
		case PixelSubdetector::PixelEndcap: {
			PXFDetId sid(detID.rawId());
			barrel = false;
			endcapType = EndcapTypes::PXF;
			layer = sid.disk() + endcapType - 1;
			side = sid.side() == 1 ? -1 : 1;
			stereo = true;
			break;
		}
		case StripSubdetector::TIB: {
			TIBDetId sid(detID.rawId());
			barrel = true;
			barrelType = BarrelTypes::TIB;
			layer = sid.layer() + barrelType - 1;
			stereo = sid.isStereo();
			break;
		}
		case StripSubdetector::TID: {
			TIDDetId sid(detID.rawId());
			barrel = false;
			endcapType = EndcapTypes::TID;
			layer = sid.wheel() + endcapType - 1;
			side = sid.side() == 1 ? -1 : 1;
			stereo = sid.isStereo();
			break;
		}
		case StripSubdetector::TOB: {
			TOBDetId sid(detID.rawId());
			barrel = true;
			barrelType = BarrelTypes::TOB;
			layer = sid.layer() + barrelType - 1;
			stereo = sid.isStereo();
			break;
		}
		case StripSubdetector::TEC: {
			TECDetId sid(detID.rawId());
			barrel = false;
			endcapType = EndcapTypes::TEC;
			layer = sid.wheel() + endcapType - 1;
			side = sid.side() == 1 ? -1 : 1;
			stereo = sid.isStereo();
			break;
		}
		default: {
			std::cout << "Error: Bad SubDetID" << std::endl;
			break;
		}
		}
	}

	LayerClassification(bool iBarrel, int iLayerSigned) {
		barrel = iBarrel;
		layer = abs(iLayerSigned);
		side = barrel ? 0 : (iLayerSigned < 0 ? -1 : 1);
		stereo = false;

		fillTypes();
	}

	LayerClassification(bool iBarrel, int iLayerUNSigned, int iSide) {
		barrel = iBarrel;
		layer = iLayerUNSigned;
		side = barrel ? 0 : iSide;
		stereo = false;

		fillTypes();
	}

	int getLayer() const {
		return (barrel ? 1 : side) * getLayerUnsigned();
	}

	int getLayerUnsigned() const {
		return layer;
	}

	int getLayerInSubDet() const {
		return layer - (barrel ? (int) barrelType : (int) endcapType) + 1;
	}

	int getSide() const {
		return side;
	}

	bool isBarrel() const {
		return barrel;
	}

	int getLadder() const {
		return ladder;
	}

	bool isStereo() const {
		return stereo;
	}

	BarrelTypes getBarrelType() const {
		return barrelType;
	}

	EndcapTypes getEndcapType() const {
		return endcapType;
	}

	//comparison operators
	inline bool operator==(const LayerClassification& rhs) const {
		return (isBarrel() == rhs.isBarrel()) && (getLayer() == rhs.getLayer());
	}

	inline bool operator!=(const LayerClassification& rhs) const {
		return !operator==(rhs);
	}

	inline bool operator<(const LayerClassification& rhs) const {
		if (isBarrel() == rhs.isBarrel()) {
			if (getSide() == rhs.getSide())
				return getLayerUnsigned() < rhs.getLayerUnsigned();
			else
				return ((int) getEndcapType()) < ((int) rhs.getEndcapType());
		} else {
			if (isBarrel())
				return ((int) getBarrelType()) < ((int) rhs.getEndcapType());
			else
				return ((int) getEndcapType())
						< (((int) rhs.getBarrelType()) - 2);
		}
	}

	inline bool operator>(const LayerClassification& rhs) const {
		return !operator<(rhs) && !operator==(rhs);
	}
	inline bool operator<=(const LayerClassification& rhs) const {
		return !operator>(rhs);
	}
	inline bool operator>=(const LayerClassification& rhs) const {
		return !operator<(rhs);
	}

	//arithmetic operators
	void operator++() {
		++layer;
		fillTypes();
	}

	LayerClassification operator+(int n) const {
		return LayerClassification(barrel, layer + n, side);
	}

	tLayers getCompatible() const {
		tLayers result;

		LayerClassification nextLayer = *this + 1;
		if ((nextLayer.isBarrel()
				&& nextLayer.barrelType != BarrelTypes::bUNKNOWN)
				|| (!nextLayer.isBarrel()
						&& nextLayer.endcapType != EndcapTypes::eUNKNOWN))
			result.push_back(nextLayer);

		//TODO define barrel endcap transistions

		return result;
	}

	bool isConsecutive(const LayerClassification& o) const {

		if (getSide() == o.getSide()) //both in either in barrel or forward or backward
			return std::abs(getLayerUnsigned() - o.getLayerUnsigned()) <= 1; //only one layer difference

		//TODO define valid barrel endcap transistions

		return false;

	}

	static const std::vector<LayerClassification> getBarrelLayers() {
		std::vector<LayerClassification> result;

		for (LayerClassification layer(true,
				LayerClassification::BarrelTypes::PXB);
				layer.getLayerUnsigned()
						< LayerClassification::BarrelTypes::bUNKNOWN; ++layer) {
			result.push_back(layer);
		}

		return result;
	}

	static const std::vector<LayerClassification> getForwardLayers() {
		std::vector<LayerClassification> result;

		for (LayerClassification layer(false,
				LayerClassification::EndcapTypes::PXF, +1);
				layer.getLayerUnsigned()
						< LayerClassification::EndcapTypes::eUNKNOWN; ++layer) {
			result.push_back(layer);
		}

		return result;
	}

	static const std::vector<LayerClassification> getBackwardLayers() {
		std::vector<LayerClassification> result;

		for (LayerClassification layer(false,
				LayerClassification::EndcapTypes::PXF, -1);
				layer.getLayerUnsigned()
						< LayerClassification::EndcapTypes::eUNKNOWN; ++layer) {
			result.push_back(layer);
		}

		return result;
	}

private:
	bool barrel; // is detector in barrel or endcap
	int layer; // for barrel: layer; for endcap: disk or wheel
	int side; // barrel = 0; endcap: pos or neg z side
	int ladder; // for pixel elements
	bool stereo; //is stereo module
	BarrelTypes barrelType;
	EndcapTypes endcapType;

	void fillTypes() {
		if (barrel) {
			endcapType = EndcapTypes::eUNKNOWN;
			if (layer < BarrelTypes::TIB)
				barrelType = BarrelTypes::PXB;
			else if (layer < BarrelTypes::TOB)
				barrelType = BarrelTypes::TIB;
			else if (layer < BarrelTypes::bUNKNOWN)
				barrelType = BarrelTypes::TOB;
			else
				barrelType = BarrelTypes::bUNKNOWN;
		} else {
			barrelType = BarrelTypes::bUNKNOWN;
			if (layer < EndcapTypes::TID)
				endcapType = EndcapTypes::PXF;
			else if (layer < EndcapTypes::TEC)
				endcapType = EndcapTypes::TID;
			else if (layer < EndcapTypes::eUNKNOWN)
				endcapType = EndcapTypes::TEC;
			else
				endcapType = EndcapTypes::eUNKNOWN;
		}
	}
};

//stream operator
std::ostream& operator<<(std::ostream& s, const LayerClassification& det);
