IMpreload('Welcome_files', 'shapeimage_21', '0');
var MINIMUM_FONT = "10";
var UNITS = "";

function elementFontSize(element)
{
    var fontSize = MINIMUM_FONT; 

    if (document.defaultView)
    {
        var computedStyle = document.defaultView.getComputedStyle(element, null);
        if (computedStyle)
        {
            fontSize = computedStyle.getPropertyValue("font-size");
        }
    }
    else if (element.currentStyle)
    {
        fontSize = element.currentStyle.fontSize;
    }

    if ((UNITS.length == 0) && (fontSize != MINIMUM_FONT))
    {
        UNITS = fontSize.substring(fontSize.length - 2, fontSize.length)
    }

    return parseFloat(fontSize);
}

function adjustFontSizeIfTooBig(idOfElement)
{
    var oTextBoxOuterDiv;
    var oTextBoxMiddleDiv;
    var oTextBoxInnerDiv;
    var oTextBoxOuterDiv = document.getElementById(idOfElement);
    
    if (oTextBoxOuterDiv)
    {
        oTextBoxMiddleDiv = getChildOfType(oTextBoxOuterDiv, "DIV", 0);
        if (oTextBoxMiddleDiv)
        {
            oTextBoxInnerDiv = getChildOfType(oTextBoxMiddleDiv, "DIV", 0);
            if (oTextBoxInnerDiv)
            {
                var offsetHeight = oTextBoxInnerDiv.offsetHeight;
                var specifiedHeight = offsetHeight;
                if (oTextBoxMiddleDiv.style.height != "")
                {
                    specifiedHeight = parseFloat(oTextBoxMiddleDiv.style.height);
                }
                else if (oTextBoxOuterDiv.style.height != "")
                {
                    specifiedHeight = parseFloat(oTextBoxOuterDiv.style.height);
                }
                if (offsetHeight > specifiedHeight)
                {
                    var smallestFontSize = 200;
                    
                    var aParaChildren = getParaDescendants(oTextBoxInnerDiv);
                    var oneLine = false;
                    for (i = 0; i < aParaChildren.length; i++)
                    {
                        var oParagraphDiv = aParaChildren[i];
                        var lineHeight = elementLineHeight(oParagraphDiv);
                        oneLine = oneLine || (lineHeight * 1.5 >= specifiedHeight);
                        if (oParagraphDiv.nodeName == "DIV")
                        {
                            var fontSize = elementFontSize(oParagraphDiv);
                            smallestFontSize = Math.min( smallestFontSize, fontSize );
                            for (j = 0; j < oParagraphDiv.childNodes.length; j++)
                            {
                                var oSpan = oParagraphDiv.childNodes[j];
                                if ((oSpan.nodeName == "SPAN") || (oSpan.nodeName == "A"))
                                {
                                    fontSize = elementFontSize(oSpan);
                                    smallestFontSize = Math.min( smallestFontSize, fontSize );
                                }
                            }
                        }
                    }
                    var minimum = parseFloat(MINIMUM_FONT);
                    
                    var count = 0
                    while ((smallestFontSize > minimum) && (offsetHeight > specifiedHeight) && (count < 10))
                    {
                        ++ count;
                        if (oneLine)
                        {
                            var oldWidth = parseInt(oTextBoxOuterDiv.style.width);
                            oTextBoxInnerDiv.style.width =
                                "" + oldWidth * Math.pow(1.05, count) + "px";
                        }
                        else
                        {
                            var scale = Math.max(0.95, minimum / smallestFontSize);
                            
                            for (i = 0; i < aParaChildren.length; i++)
                            {
                                var oParagraphDiv = aParaChildren[i];
                                if (oParagraphDiv.nodeName == "DIV")
                                {
                                    var paraFontSize = elementFontSize(oParagraphDiv) * scale;
                                    var paraLineHeight = elementLineHeight(oParagraphDiv) * scale;
                                    for (j = 0; j < oParagraphDiv.childNodes.length; j++)
                                    {
                                        var oSpan = oParagraphDiv.childNodes[j];
                                        if ((oSpan.nodeName == "SPAN") || (oSpan.nodeName == "A"))
                                        {
                                            var spanFontSize = elementFontSize(oSpan) * scale;
                                            var spanLineHeight = elementLineHeight(oSpan) * scale;
                                            oSpan.style.fontSize = spanFontSize + UNITS;
                                            oSpan.style.lineHeight = spanLineHeight + UNITS;
                                            smallestFontSize = Math.min( smallestFontSize, spanFontSize );
                                        }
                                    }
                                    oParagraphDiv.style.fontSize = paraFontSize + UNITS;
                                    oParagraphDiv.style.lineHeight = paraLineHeight + UNITS;
                                    smallestFontSize = Math.min( smallestFontSize, paraFontSize );
                                }
                            }
                        }
                        
                        offsetHeight = oTextBoxInnerDiv.offsetHeight;
                    }
                }
            }
        }
    }
}


function elementLineHeight(element)
{
    var lineHeight = MINIMUM_FONT; 
    
    if (document.defaultView)
    {
        var computedStyle = document.defaultView.getComputedStyle(element, null);
        if (computedStyle)
        {
            lineHeight = computedStyle.getPropertyValue("line-height");
        }
    }
    else if (element.currentStyle)
    {
        lineHeight = element.currentStyle.lineHeight;
    }
    
    if ((UNITS.length == 0) && (lineHeight != MINIMUM_FONT))
    {
        UNITS = lineHeight.substring(lineHeight.length - 2, lineHeight.length)
    }
    
    return parseFloat(lineHeight);
}

function adjustLineHeightIfTooBig(idOfElement)
{
    var oTextBoxOuterDiv;
    var oTextBoxMiddleDiv;
    var oTextBoxInnerDiv;
    var oTextBoxOuterDiv = document.getElementById(idOfElement);
    
    if (oTextBoxOuterDiv)
    {
        oTextBoxMiddleDiv = getChildOfType(oTextBoxOuterDiv, "DIV", 0);
        if (oTextBoxMiddleDiv)
        {
            oTextBoxInnerDiv = getChildOfType(oTextBoxMiddleDiv, "DIV", 0);
            if (oTextBoxInnerDiv)
            {
                var offsetHeight = oTextBoxInnerDiv.offsetHeight;
                var specifiedHeight = offsetHeight;
                if (oTextBoxMiddleDiv.style.height != "")
                {
                    specifiedHeight = parseFloat(oTextBoxMiddleDiv.style.height);
                }
                else if (oTextBoxOuterDiv.style.height != "")
                {
                    specifiedHeight = parseFloat(oTextBoxOuterDiv.style.height);
                }
                if (offsetHeight > specifiedHeight)
                {
                    var adjusted = true;
                    var count = 0;
                    while ((adjusted) && (offsetHeight > specifiedHeight) && (count < 10))
                    {
                        adjusted = false;
                        ++ count;
                        
                        var aParaChildren = getParaDescendants(oTextBoxInnerDiv);
                        for (i = 0; i < aParaChildren.length; i++)
                        {
                            var oParagraphDiv = aParaChildren[i];
                            if (oParagraphDiv.nodeName == "DIV")
                            {
                                var fontSize = elementFontSize(oParagraphDiv);
                                var lineHeight = elementLineHeight(oParagraphDiv) * 0.95;
                                if (lineHeight >= (fontSize * 1.1))
                                {
                                    oParagraphDiv.style.lineHeight = lineHeight + UNITS;
                                    adjusted = true;
                                }
                                
                                
                                
                                for (j = 0; j < oParagraphDiv.childNodes.length; j++)
                                {
                                    var oSpan = oParagraphDiv.childNodes[j];
                                    if ((oSpan.nodeName == "SPAN") || (oSpan.nodeName == "A"))
                                    {
                                        var fontSize = elementFontSize(oSpan);
                                        var lineHeight = elementLineHeight(oSpan) * 0.95;
                                        if (lineHeight >= (fontSize * 1.1))
                                        {
                                            oSpan.style.lineHeight = lineHeight + UNITS;
                                            var adjusted = true;
                                        }
                                    }
                                }
                            }
                        }
                        
                        offsetHeight = oTextBoxInnerDiv.offsetHeight;
                    }
                }
            }
        }
    }
}

var smallTransparentGif = "";
function fixupIEPNG(strImageID, transparentGif) 
{
    smallTransparentGif = transparentGif;
    if (windowsInternetExplorer && (browserVersion < 7))
    {
        var img = document.getElementById(strImageID);
        if (img)
        {
            var src = img.src;
            img.style.filter = "progid:DXImageTransform.Microsoft.AlphaImageLoader(src='" + src + "', sizingMethod='scale')";
            img.src = transparentGif;
            img.attachEvent("onpropertychange", imgPropertyChanged);
        }
    }
}

function fixupIEPNGBG(oBlock) 
{
    if (oBlock)
    {
        var currentBGImage = oBlock.currentStyle.backgroundImage;
        var currentBGRepeat = oBlock.currentStyle.backgroundRepeat;
        var urlStart = currentBGImage.indexOf('url(');
        var urlEnd = currentBGImage.indexOf(')', urlStart);
        var imageURL = currentBGImage.substring(urlStart + 4, urlEnd);

        if (imageURL.charAt(0) == '"')
        {
            imageURL = imageURL.substring(1);
        }
        
        if (imageURL.charAt(imageURL.length - 1) == '"')
        {
            imageURL = imageURL.substring(0, imageURL.length - 1);
        }

        var overrideRepeat = false;

        var filterStyle =
            "progid:DXImageTransform.Microsoft.AlphaImageLoader(src='" +
            imageURL +
            "', sizingMethod='crop');";

        if (RegExp("/C[0-9A-F]{8}.png$").exec(imageURL) != null)
        {
            filterStyle =
                "progid:DXImageTransform.Microsoft.AlphaImageLoader(src='" +
                imageURL +
                "', sizingMethod='scale');";

            overrideRepeat = true;
        }

        var backgroundImage = new Image();
        backgroundImage.src = imageURL;
        var tileWidth = backgroundImage.width;
        var tileHeight = backgroundImage.height; 
        
        var blockWidth = 0;
        var blockHeight = 0;
        if (oBlock.style.width)
        {
            blockWidth = parseInt(oBlock.style.width);
        }
        else
        {
            blockWidth = oBlock.offsetWidth;
        }
        if (oBlock.style.height)
        {
            blockHeight = parseInt(oBlock.style.height);
        }
        else
        {
            blockHeight = oBlock.offsetHeight;
        }

        if ((blockWidth == 0) || (blockHeight == 0))
        {
            return;
        }
        
        var wholeRows = 1;
        var wholeCols = 1;
        var extraHeight = 0;
        var extraWidth = 0;
        
        if ((currentBGRepeat.indexOf("no-repeat") != -1) ||
              ((tileWidth == 0) && (tileHeight == 0)) ||
              overrideRepeat)
        {
            tileWidth = blockWidth;
            tileHeight = blockHeight;

        }
        else if ((currentBGRepeat.indexOf("repeat-x") != -1) ||
              (tileHeight == 0))
        {
            wholeCols = Math.floor(blockWidth / tileWidth);
            extraWidth = blockWidth - (tileWidth * wholeCols);
            tileHeight = blockHeight;

        }
        else if (currentBGRepeat.indexOf("repeat-y") != -1)
        {
            wholeRows = Math.floor(blockHeight / tileHeight);
            extraHeight = blockHeight - (tileHeight * wholeRows);
            tileWidth = blockWidth;

        }
        else
        {
            wholeCols = Math.floor(blockWidth / tileWidth);
            wholeRows = Math.floor(blockHeight / tileHeight);
            extraWidth = blockWidth - (tileWidth * wholeCols);
            extraHeight = blockHeight - (tileHeight * wholeRows);
        }
        
        var wrappedContent = document.createElement("div");
        wrappedContent.style.position = "relative";
        wrappedContent.style.zIndex = "1";
        wrappedContent.style.left = "0px";
        wrappedContent.style.top = "0px";
        if (!isNaN(parseInt(oBlock.style.width)))
        {
            wrappedContent.style.width = "" + blockWidth + "px";
        }
        if (!isNaN(parseInt(oBlock.style.height)))
        {
            wrappedContent.style.height = "" + blockHeight + "px";
        }
        var pngBGFixIsWrappedContentEmpty = true;
        while (oBlock.hasChildNodes())
        {
            if (oBlock.firstChild.nodeType == 3)
            {
                if (RegExp("^ *$").exec(oBlock.firstChild.data) == null)
                {
                    pngBGFixIsWrappedContentEmpty = false;
                }
            }
            else
            {
                pngBGFixIsWrappedContentEmpty = false;
            }
            wrappedContent.appendChild(oBlock.firstChild);
        }
        if (pngBGFixIsWrappedContentEmpty)
        {
            wrappedContent.style.lineHeight = "0px";
        }
        
        var newMarkup = "";
        for (var currentRow = 0; 
             currentRow < wholeRows; 
             currentRow++)
        {
            for (currentCol = 0; 
                 currentCol < wholeCols; 
                 currentCol++)
            {
                newMarkup += "<div style=" +
                        "\"position: absolute; line-height: 0px; " +
                        "width: " + tileWidth + "px; " +
                        "height: " + tileHeight + "px; " +
                        "left:" + currentCol *  tileWidth + "px; " +
                        "top:" + currentRow *  tileHeight + "px; " +
                        "filter:" + filterStyle + 
                        "\" > </div>";
            }
            
            if (extraWidth != 0)
            {
                newMarkup += "<div style=" +
                        "\"position: absolute; line-height: 0px; " +
                        "width: " + extraWidth + "px; " +
                        "height: " + tileHeight + "px; " +
                        "left:" + currentCol *  tileWidth + "px; " +
                        "top:" + currentRow *  tileHeight + "px; " +
                        "filter:" + filterStyle + 
                        "\" > </div>";
            }
        }
        
        if (extraHeight != 0)
        {
            for (currentCol = 0; 
                 currentCol < wholeCols; 
                 currentCol++)
            {
                newMarkup += "<div style=" +
                        "\"position: absolute; line-height: 0px; " +
                        "width: " + tileWidth + "px; " +
                        "height: " + extraHeight + "px; " +
                        "left:" + currentCol *  tileWidth + "px; " +
                        "top:" + currentRow *  tileHeight + "px; " +
                        "filter:" + filterStyle + 
                        "\" > </div>";
            }
            
            if (extraWidth != 0)
            {
                newMarkup += "<div style=" +
                        "\"position: absolute; line-height: 0px; " +
                        "width: " + extraWidth + "px; " +
                        "height: " + extraHeight + "px; " +
                        "left:" + currentCol *  tileWidth + "px; " +
                        "top:" + currentRow *  tileHeight + "px; " +
                        "filter:" + filterStyle + 
                        "\" > </div>";
            }
        }
        oBlock.innerHTML = newMarkup;

        oBlock.appendChild(wrappedContent);
        oBlock.style.background= "";
    }
}

function fixupAllIEPNGBGs()
{
    if (windowsInternetExplorer && (browserVersion < 7))
    {
        try
        {
            var oDivNodes = document.getElementsByTagName('DIV');
            for (var iIndex=0; iIndex<oDivNodes.length; iIndex++)
            {
                var oNode = oDivNodes.item(iIndex);
                if (oNode.currentStyle &&
                    oNode.currentStyle.backgroundImage &&
                    (oNode.currentStyle.backgroundImage.indexOf('url(') != -1) &&
                    (oNode.currentStyle.backgroundImage.indexOf('.png")') != -1))
                {
                    fixupIEPNGBG(oNode);
                }
            }
        }
        catch (e)
        {
        }
    }
}

function getChildOfType(oParent, sNodeName, requestedIndex)
{
    var childrenOfType = oParent.getElementsByTagName(sNodeName);
    return (requestedIndex < childrenOfType.length) ?
           childrenOfType.item(requestedIndex) : null;
}

function onPageLoad()
{
    detectBrowser();
    adjustLineHeightIfTooBig("id2");
    adjustFontSizeIfTooBig("id2");
    adjustLineHeightIfTooBig("id6");
    adjustFontSizeIfTooBig("id6");
    adjustLineHeightIfTooBig("id8");
    adjustFontSizeIfTooBig("id8");
    adjustLineHeightIfTooBig("id10");
    adjustFontSizeIfTooBig("id10");
    adjustLineHeightIfTooBig("id12");
    adjustFontSizeIfTooBig("id12");
    adjustLineHeightIfTooBig("id14");
    adjustFontSizeIfTooBig("id14");
    adjustLineHeightIfTooBig("id16");
    adjustFontSizeIfTooBig("id16");
    adjustLineHeightIfTooBig("id18");
    adjustFontSizeIfTooBig("id18");
    adjustLineHeightIfTooBig("id31");
    adjustFontSizeIfTooBig("id31");
    adjustLineHeightIfTooBig("id32");
    adjustFontSizeIfTooBig("id32");
    adjustLineHeightIfTooBig("id33");
    adjustFontSizeIfTooBig("id33");
    adjustLineHeightIfTooBig("id36");
    adjustFontSizeIfTooBig("id36");
    adjustLineHeightIfTooBig("id38");
    adjustFontSizeIfTooBig("id38");
    adjustLineHeightIfTooBig("id39");
    adjustFontSizeIfTooBig("id39");
    adjustLineHeightIfTooBig("id41");
    adjustFontSizeIfTooBig("id41");
    adjustLineHeightIfTooBig("id42");
    adjustFontSizeIfTooBig("id42");
    adjustLineHeightIfTooBig("id43");
    adjustFontSizeIfTooBig("id43");
    adjustLineHeightIfTooBig("id44");
    adjustFontSizeIfTooBig("id44");
    adjustLineHeightIfTooBig("id45");
    adjustFontSizeIfTooBig("id45");
    adjustLineHeightIfTooBig("id46");
    adjustFontSizeIfTooBig("id46");
    adjustLineHeightIfTooBig("id47");
    adjustFontSizeIfTooBig("id47");
    adjustLineHeightIfTooBig("id48");
    adjustFontSizeIfTooBig("id48");
    adjustLineHeightIfTooBig("id65");
    adjustFontSizeIfTooBig("id65");
    adjustLineHeightIfTooBig("id68");
    adjustFontSizeIfTooBig("id68");
    adjustLineHeightIfTooBig("id69");
    adjustFontSizeIfTooBig("id69");
    adjustLineHeightIfTooBig("id70");
    adjustFontSizeIfTooBig("id70");
    adjustLineHeightIfTooBig("id93");
    adjustFontSizeIfTooBig("id93");
    adjustLineHeightIfTooBig("id94");
    adjustFontSizeIfTooBig("id94");
    adjustLineHeightIfTooBig("id95");
    adjustFontSizeIfTooBig("id95");
    fixupAllIEPNGBGs();
    fixupIEPNG("id1", "Welcome_files/transparent.gif");
    fixupIEPNG("id3", "Welcome_files/transparent.gif");
    fixupIEPNG("id4", "Welcome_files/transparent.gif");
    fixupIEPNG("id5", "Welcome_files/transparent.gif");
    fixupIEPNG("id7", "Welcome_files/transparent.gif");
    fixupIEPNG("id9", "Welcome_files/transparent.gif");
    fixupIEPNG("id11", "Welcome_files/transparent.gif");
    fixupIEPNG("id13", "Welcome_files/transparent.gif");
    fixupIEPNG("id15", "Welcome_files/transparent.gif");
    fixupIEPNG("id17", "Welcome_files/transparent.gif");
    fixupIEPNG("id19", "Welcome_files/transparent.gif");
    fixupIEPNG("id20", "Welcome_files/transparent.gif");
    fixupIEPNG("id21", "Welcome_files/transparent.gif");
    fixupIEPNG("id22", "Welcome_files/transparent.gif");
    fixupIEPNG("id23", "Welcome_files/transparent.gif");
    fixupIEPNG("id24", "Welcome_files/transparent.gif");
    fixupIEPNG("id25", "Welcome_files/transparent.gif");
    fixupIEPNG("id26", "Welcome_files/transparent.gif");
    fixupIEPNG("id27", "Welcome_files/transparent.gif");
    fixupIEPNG("id28", "Welcome_files/transparent.gif");
    fixupIEPNG("id29", "Welcome_files/transparent.gif");
    fixupIEPNG("id30", "Welcome_files/transparent.gif");
    fixupIEPNG("shapeimage_21_link_0", "Welcome_files/transparent.gif");
    fixupIEPNG("shapeimage_21", "Welcome_files/transparent.gif");
    fixupIEPNG("id34", "Welcome_files/transparent.gif");
    fixupIEPNG("id35", "Welcome_files/transparent.gif");
    fixupIEPNG("id37", "Welcome_files/transparent.gif");
    fixupIEPNG("id40", "Welcome_files/transparent.gif");
    fixupIEPNG("id49", "Welcome_files/transparent.gif");
    fixupIEPNG("id50", "Welcome_files/transparent.gif");
    fixupIEPNG("id51", "Welcome_files/transparent.gif");
    fixupIEPNG("id52", "Welcome_files/transparent.gif");
    fixupIEPNG("id53", "Welcome_files/transparent.gif");
    fixupIEPNG("id54", "Welcome_files/transparent.gif");
    fixupIEPNG("id55", "Welcome_files/transparent.gif");
    fixupIEPNG("id56", "Welcome_files/transparent.gif");
    fixupIEPNG("id57", "Welcome_files/transparent.gif");
    fixupIEPNG("id58", "Welcome_files/transparent.gif");
    fixupIEPNG("id59", "Welcome_files/transparent.gif");
    fixupIEPNG("id60", "Welcome_files/transparent.gif");
    fixupIEPNG("id61", "Welcome_files/transparent.gif");
    fixupIEPNG("id62", "Welcome_files/transparent.gif");
    fixupIEPNG("id63", "Welcome_files/transparent.gif");
    fixupIEPNG("id64", "Welcome_files/transparent.gif");
    fixupIEPNG("id66", "Welcome_files/transparent.gif");
    fixupIEPNG("id67", "Welcome_files/transparent.gif");
    fixupIEPNG("id71", "Welcome_files/transparent.gif");
    fixupIEPNG("id72", "Welcome_files/transparent.gif");
    fixupIEPNG("id73", "Welcome_files/transparent.gif");
    fixupIEPNG("id74", "Welcome_files/transparent.gif");
    fixupIEPNG("id75", "Welcome_files/transparent.gif");
    fixupIEPNG("id76", "Welcome_files/transparent.gif");
    fixupIEPNG("id77", "Welcome_files/transparent.gif");
    fixupIEPNG("id78", "Welcome_files/transparent.gif");
    fixupIEPNG("id79", "Welcome_files/transparent.gif");
    fixupIEPNG("id80", "Welcome_files/transparent.gif");
    fixupIEPNG("id81", "Welcome_files/transparent.gif");
    fixupIEPNG("id82", "Welcome_files/transparent.gif");
    fixupIEPNG("id83", "Welcome_files/transparent.gif");
    fixupIEPNG("id84", "Welcome_files/transparent.gif");
    fixupIEPNG("id85", "Welcome_files/transparent.gif");
    fixupIEPNG("id86", "Welcome_files/transparent.gif");
    fixupIEPNG("id87", "Welcome_files/transparent.gif");
    fixupIEPNG("id88", "Welcome_files/transparent.gif");
    fixupIEPNG("id89", "Welcome_files/transparent.gif");
    fixupIEPNG("id90", "Welcome_files/transparent.gif");
    fixupIEPNG("id91", "Welcome_files/transparent.gif");
    fixupIEPNG("id92", "Welcome_files/transparent.gif");
    fixupIEPNG("id96", "Welcome_files/transparent.gif");
    fixupIEPNG("id97", "Welcome_files/transparent.gif");
    fixupIEPNG("id98", "Welcome_files/transparent.gif");
    fixupIEPNG("id99", "Welcome_files/transparent.gif");
    return true;
}

function getParaDescendants(oAncestor)
{
    var oParaDescendants = new Array();
    var oPotentialParagraphs = oAncestor.getElementsByTagName('DIV');
    for (var iIndex=0; iIndex<oPotentialParagraphs.length; iIndex++)
    {
        var oNode = oPotentialParagraphs.item(iIndex);
        if (oNode.className.lastIndexOf('paragraph') != -1)
        {
            oParaDescendants.push(oNode);
        }
    }
    return oParaDescendants;
}

function NBmouseover(index)
{
    var normal = document.getElementById("navbar_"+index+"_normal");
    var rollover = document.getElementById("navbar_"+index+"_rollover");
    if (normal && rollover)
    {
        normal.style.visibility = "hidden";
        rollover.style.visibility = "visible";
    }
    return true;
}

function NBmouseout(index)
{
    var normal = document.getElementById("navbar_"+index+"_normal");
    var rollover = document.getElementById("navbar_"+index+"_rollover");
    if (normal && rollover)
    {
        normal.style.visibility = "visible";
        rollover.style.visibility = "hidden";
    }
    return true;
}

function IMpreload(path, name, areaIndex)
{
    var rolloverName = name+'_rollover_'+areaIndex;
    var rolloverPath = path+'/'+rolloverName+'.png';
    self[rolloverName] = new Image();
    self[rolloverName].src = rolloverPath;

    var linkName = name+'_link_'+areaIndex;
    var linkPath = path+'/'+linkName+'.png';
    self[linkName] = new Image();
    self[linkName].src = linkPath;
    return true;
}

function IMmouseover(name, areaIndex)
{
    var rolloverName = name+'_rollover_'+areaIndex;
    var linkName = name+'_link_'+areaIndex;
    var img  = document.getElementById(linkName);
    if (img)
    {
        img.src = self[rolloverName].src;
    }
    return true;
}

function IMmouseout(name, areaIndex)
{
    var linkName = name+'_link_'+areaIndex;
    var img  = document.getElementById(linkName);
    if (img)
    {
        img.src = self[linkName].src;
    }
    return true;
}

var windowsInternetExplorer = false;
var browserVersion = 0;
function detectBrowser()
{
    windowsInternetExplorer = false;
    var appVersion = navigator.appVersion;
    if ((appVersion.indexOf("MSIE") != -1) &&
        (appVersion.indexOf("Macintosh") == -1))
    {
        var temp = appVersion.split("MSIE");
        browserVersion = parseFloat(temp[1]);
        windowsInternetExplorer = true;
    }
}

var inImgPropertyChanged = false;
function imgPropertyChanged()
{
    if ((window.event.propertyName == "src") && (! inImgPropertyChanged))
    {
        inImgPropertyChanged = true;
        var el = window.event.srcElement;
        if (el.src != smallTransparentGif)
        {
            el.filters.item(0).src = el.src;
            el.src = smallTransparentGif;
        }
        inImgPropertyChanged = false;
    }
}

